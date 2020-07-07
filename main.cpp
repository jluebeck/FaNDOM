//Authors: Jens Luebeck (jluebeck@ucsd.edu), Siavash Raisi (sraeisid@ucsd.edu)

#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <map>
#include <ctime>
#include <future>
#include <chrono>

#include "threadsafe_queue.h"
#include "OMIO.h"
#include "OMHelper.h"
#include "fandomAligner.h"
#include "OMFilter.h"

using namespace std;

static string version = "0.1";

//Alignment variables
const int lookback = 5;
const int min_map_lab = 10;
const int min_map_len = 25000;
float aln_padding = 1000;
bool multimap_mols = false;
bool partial_alignment = true;
int score_limit = 5000;

//Data filtering
int n_threads = 1;
int min_aln_len = 6;
float aln_prop_thresh_to_remap = 0.8;
double aln_len_thresh_to_remap = 25000.;

//bool subsect = false;
//map<int,vector<pair<int,int>>> bed_data;
ofstream logfile;

//------------------------------------------------------------------------------------------
//Run Fandom alignment on a batch of seeds
map<int, vector<Alignment>> run_aln(map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &mols,
        map<int, vector<seedData>> &mol_seed_data, bool partial_mode) {
    map<int, vector<Alignment>> results;

    //iterate over the molecules
    for (auto &item: mol_seed_data) {
        int mol_id = item.first;
        vector<seedData> curr_mol_seed_data = item.second;
        double best_score = 0;
        vector<Alignment> mol_alns;
        vector<double> mol_vect = mols[mol_id];
        //TODO: This following function in OMHelper creates a collection of alignment windows from the seeds it is given. This step can be optimized
        map<int, vector<seedData>> curr_mol_seeds = mol_seeds_to_aln_regions(curr_mol_seed_data, ref_cmaps, mol_vect, aln_padding);

        //iterate over the different reference contigs hit by the seeds
        for (auto &x: curr_mol_seeds) {
            int ref_id = x.first;
            vector<seedData> seed_list = x.second;
            vector<double> ref_vect = ref_cmaps[ref_id];

            //iterate over the individual seeds
            for (const auto &s: seed_list) {
                int a = s.ref_aln_lb;
                int b = s.ref_aln_rb;
                //initialize storage for alignment
                vector<vector<double>> S((b - a), vector<double>(mol_vect.size() - 1, 0));
                vector<vector<pair<int, int>>> previous((b - a),
                                                        vector<pair<int, int>>(mol_vect.size() - 1, {-1, -1}));

                //perform alignment
                dp_aln(S, previous, mol_vect, ref_vect, a, b, lookback);
                pair<int,int> max_pair =  get_max_pair(S);
                double curr_score = S[max_pair.first][max_pair.second];

                //check if the alignment passes cutoffs
                if (multimap_mols || curr_score > best_score || partial_mode) {
                    Alignment curr_aln = dp_backtracking(S, previous, max_pair, a, ref_id, mol_id);
                    curr_aln.seed_num = s.seed_num;
                    if (partial_mode) {
                        //partial alignment confidence check
                        if (! partial_confidence_check(curr_aln, ref_cmaps)) {
                            continue;
                        }
                        curr_aln.is_partial = 1;
                    }
                    if (curr_aln.alignment.size() < min_aln_len ||
                        curr_score / curr_aln.alignment.size() < score_limit) {
                        continue;
                    }
                    mol_alns.push_back(curr_aln);

                }
                if (curr_score > best_score) {
                    best_score = curr_score;
                }
            }
        }

        //go over alignments and check if the alignments are secondary alignments and set the status
        int is_multimapped = (int) (mol_alns.size() > 1);
        cout << mol_id << " " << lround(best_score) << "\n";
        for (auto &curr_aln: mol_alns) {
            curr_aln.is_multimapped = is_multimapped;
            if (get<2>(curr_aln.alignment.back()) < lround(best_score)) {
                curr_aln.is_secondary = 1;
            } else {
                cout << " hit " << get<2>(curr_aln.alignment.back()) << "\n";
            }
        }
        results[mol_id] = mol_alns;
    }
    return results;
}

//run the filter and call the alignment
map<int, vector<Alignment>> filt_and_aln(int thread_num, map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &mols,
                                        map<int, dis_to_index> &ref_DTI, map<int, int> &ref_lens,
                                        threadsafe_queue<int> &mol_id_queue, bool partial_mode) {
    int counter = 0;
    vector<query> qq;
    map<int, vector<seedData>> seed_batch;
    map<int, vector<Alignment>> curr;
    map<int, vector<Alignment>> result;

    //while molecules in queue pop off
    while (true) {
        int q_id = mol_id_queue.pop();
        if (!q_id) {
            break;
        } else {
            int qsize = mol_id_queue.size();
            if (qsize % 10000 == 0) {
                cout << qsize << " entries left in queue\n";
            }
            vector<double> q_posns = mols[q_id];
            query q (q_id, calculate_genomic_distance(mols[q_id]),mols[q_id].size() - 1,0);
            qq.push_back(q);
            counter = counter + q.length + w;
        }
        //if the storage for the seed list isn't full
        if (counter > 5000) {
            seed_batch = OMFilter(qq, thread_num, ref_DTI, ref_lens, ref_cmaps, mols, partial_mode);
            curr = run_aln(ref_cmaps, mols, seed_batch, partial_mode);
            //////Log
            for (auto &key : curr){
                logfile<< "used in full sort " << key.first<<"\n";
            }
            //////
            result.insert(curr.begin(), curr.end());
            counter = 0;
            qq.clear();
        }
    }
    //cleanup the case where the last element handled didn't fill up the seed list
    if (! qq.empty()) {
        seed_batch = OMFilter(qq, thread_num, ref_DTI, ref_lens, ref_cmaps, mols, partial_mode);
        curr = run_aln(ref_cmaps, mols, seed_batch, partial_mode);
        //////Log
        for (auto &key : curr){
            cout<< "used in cleanup sort " << key.first<<endl;
        }
        ////////
        result.insert(curr.begin(), curr.end());
    }
    return result;
}
//-----------------------------------------------------------------------
//check filenames
bool hasEnding (string const& fullString, string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}
//handle args
tuple<string,string,string,string,string,string> parse_args(int argc, char *argv[]) {
    string ref_cmap_file, queryfile, bedfile, keyfile, sample_name, outfmt;
    outfmt = "xmap";
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]).rfind("-t=", 0) == 0) {
            n_threads = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if (string(argv[i]).rfind("-r=", 0) == 0) {
            ref_cmap_file = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-q=", 0) == 0) {
            queryfile = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-sname=", 0) == 0) {
            sample_name = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-padding=", 0) == 0) {
            aln_padding = stof(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-w=", 0) == 0)) {
            w = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-band_width=", 0) == 0)) {
            band_width = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-tolerance=", 0) == 0)) {
            tolerance = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-rank=", 0) == 0)) {
            ranked = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-threshold=", 0) == 0)) {
            threshold = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if ((string(argv[i]).rfind("-outfmt=", 0) == 0)) {
            outfmt = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if ((string(argv[i]).rfind("-no_partial", 0) == 0)) {
//                Partial_aligning = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
            partial_alignment = false;

        } else if ((string(argv[i]).rfind("-multimap", 0) == 0)) {
            multimap_mols = true;

        } else if ((string(argv[i]).rfind("-version", 0) == 0)) {
            cout << "FaNDOM version " << version << "\n";
            exit(0);
        }
    }
    if (ref_cmap_file.empty()) {
        cerr << "Reference genome .cmap file [-r=] unspecified\n";
        exit(1);
    } else if (queryfile.empty()) {
        cerr << "Query file [-q=] unspecified\n";
        exit(1);
    } else if (sample_name.empty()) {
        cerr << "Sample output name [-s=] unspecified\n";
    }

    for (const string &fname: {ref_cmap_file,queryfile}) {
        if (!file_exists(fname)) {
            cerr << fname << " does not exist\n";
            exit(1);
        }
    }
    transform(outfmt.begin(), outfmt.end(), outfmt.begin(), ::tolower);
    if (outfmt != "xmap" && outfmt != "fda") {
        cerr << "-outfmt must be xmap or fda";
        exit(1);
    }
    return make_tuple(ref_cmap_file,queryfile,bedfile,keyfile,sample_name,outfmt);
}
/*
 *MAIN
 */
int main (int argc, char *argv[]) {
    if (argc < 2) {
        cout << "wrong number of arguments" << endl;
        exit(0);
    }
    //timing
    clock_t begin = clock();
    chrono::steady_clock::time_point beginWall = chrono::steady_clock::now();

    //-----------------------------
    //get args
    string ref_cmap_file, queryfile, bedfile, keyfile, sample_name, outfmt;
    tie(ref_cmap_file,queryfile,bedfile,keyfile, sample_name, outfmt) = parse_args(argc,argv);
    string log_dir = sample_name + "_log";

    logfile.open(log_dir);
    logfile<<"KIR"<<endl;
    cout << "reference genome cmap file: " << ref_cmap_file << "\n";
    cout << "query file: " << queryfile << "\n";
    cout << "\nNumber of threads for filtering & alignment: " << n_threads << "\n";
    //------------------------
    //parse inputs
    //parse ref genome
    chrono::steady_clock::time_point readWallS = chrono::steady_clock::now();

    map<int,vector<double>> ref_genome_unrev = parse_cmap(ref_cmap_file);
    cout << "Finished reading the ref\n";
    map<int,vector<double>> ref_cmaps = make_reverse_cmap(ref_genome_unrev);

    cout << "Reading queries\n";
    map<int, vector<double>> mol_maps;
    if (hasEnding(queryfile,".bnx")) {
        mol_maps = parse_bnx(queryfile);
    } else if (hasEnding(queryfile,".cmap")) {
        mol_maps = parse_cmap(queryfile);
    } else {
        perror("query file had unrecognized file extension\n");
        return 1;
    }

    cout << "Parsed query. Preprocessing...\n";
    chrono::steady_clock::time_point readWallE = chrono::steady_clock::now();

    //-----------------------------------------------------------
    //Preprocessing data for alignments

    //map of number of labels per entry
    map<int, int> ref_num_to_length = calculate_length(ref_cmaps);

    //create map of genomic distance to label number
    map<int, dis_to_index> ref_DTI = genomic_distance(ref_genome_unrev);

    //make the mol_id_queue
    chrono::steady_clock::time_point ppWallS = chrono::steady_clock::now();
    filter_mols(mol_maps,min_map_lab,min_map_len);
    cout << "Got " << mol_maps.size() << " molecules passing initial filters\n";

    threadsafe_queue<int> mol_id_queue;
    for (const auto &i: mol_maps) {
        int mol_id = i.first;
        mol_id_queue.push(mol_id);
    }

    //set outfile names
    string outname = sample_name;
    cout << "Writing alignments\n";
    string partialname = outname + "_partial";
    string argstring;
    for (int i = 0; i < argc; ++i) {
        argstring += " ";
        argstring += argv[i];
    }

    ofstream full_outfile;
    ofstream partial_outfile;
    if (outfmt == "xmap") {
        outname += ".xmap";
        partialname += ".xmap";
        full_outfile.open(outname);
        write_xmap_header(full_outfile, argstring);
        if (partial_alignment) {
            partial_outfile.open(partialname);
            write_xmap_header(partial_outfile, argstring);
        }

    } else {
        outname+=".fda";
        partialname+=".fda";
        full_outfile.open(outname);
        write_fda_header(full_outfile);
        if (partial_alignment) {
            partial_outfile.open(partialname);
            write_fda_header(partial_outfile);
        }
    }
    chrono::steady_clock::time_point ppWallE = chrono::steady_clock::now();

    //------------------------------------------------------
    //launch alignments
    chrono::steady_clock::time_point alnWallS = chrono::steady_clock::now();
    cout << "Performing alignments\n";
    vector<future< map<int, vector<Alignment>>>> futs;
    vector<promise< map<int, vector<Alignment>>>> promises(n_threads);
   for (int i = 0; i < n_threads; i++) {
       futs.push_back(async(launch:: async, filt_and_aln, i, ref(ref_cmaps), ref(mol_maps), ref(ref_DTI),
               ref(ref_num_to_length), ref(mol_id_queue), false));
   }
    //-------------------------------------------------------
    //gather results FULL
    vector<Alignment> combined_results;
    vector<Alignment> current_result;
    size_t total_alns = 0;
    for (auto &f: futs) {
       map<int, vector<Alignment>> curr_result = f.get();
       for (const auto &x: curr_result) {
           current_result = x.second;
           combined_results.reserve(combined_results.size() + distance(x.second.begin(),x.second.end()));
           combined_results.insert(combined_results.end(),x.second.begin(),x.second.end());
           //write it
           if (outfmt == "xmap") {
               write_xmap_alignment(full_outfile, current_result, ref_cmaps, mol_maps, multimap_mols, total_alns+1);
           } else {
               write_fda_alignment(full_outfile, current_result, ref_cmaps, mol_maps, multimap_mols, total_alns+1);
           }
           total_alns+=x.second.size();
       }
    }

    full_outfile << flush;
    full_outfile.close();
    cout << "Finished non-partial molecule alignment. \n" << total_alns << " total alignments\n";

    //------------------------------------------------------
    //gather results PARTIAL
//    vector<Alignment> combined_results_partial;
    if(partial_alignment){
        unordered_set<int> mols_to_remap = get_remap_mol_ids(combined_results, mol_maps, aln_prop_thresh_to_remap,
                                                             aln_len_thresh_to_remap);

        cout << mols_to_remap.size() << " molecules will undergo partial-seeding.\n";
        logfile << mols_to_remap.size() << " molecules will undergo partial-seeding.\n";
        ////////////////////
        score_limit = 1000;
        dist_scale_par = 1.4;
        penalty_par = 7500;
        //Here again make thread and run partial alignments for remaining molecules
        if (!mols_to_remap.empty()) {
            cout << "Performing partial alignments\n";
            logfile<< "Performing partial alignments\n";
            threadsafe_queue<int> mol_id_queue_partial;
            for (auto &i: mols_to_remap) {
                mol_id_queue_partial.push(i);
            }
            vector<future< map<int, vector<Alignment>>>> futs_partial;
            vector<promise< map<int, vector<Alignment>>>> promises_partial(n_threads);
            for (int i = 0; i < n_threads; i++) {
                futs_partial.push_back(async(launch::async, filt_and_aln, i, ref(ref_cmaps), ref(mol_maps), ref(ref_DTI),
                                             ref(ref_num_to_length), ref(mol_id_queue_partial), true));
            }
            total_alns = 0;
            for (auto &f: futs_partial) {
                map<int, vector<Alignment>> curr_result = f.get();
                for (const auto &x: curr_result) {
                    //write it
                    current_result = x.second;
                    if (outfmt == "xmap") {
                        write_xmap_alignment(partial_outfile, current_result, ref_cmaps, mol_maps, multimap_mols, total_alns+1);
                    } else {
                        write_fda_alignment(partial_outfile, current_result, ref_cmaps, mol_maps, multimap_mols, total_alns+1);
                    }
                    total_alns+=x.second.size();
                }

            }
            cout << "Finished partial molecule alignment. \n" << total_alns << " total alignments\n";
        }
    }
    chrono::steady_clock::time_point alnWallE = chrono::steady_clock::now();

    ////////////////////
    //write output
//    chrono::steady_clock::time_point outWallS = chrono::steady_clock::now();
//    cout << "Writing fulls\n";
//    write_xmap_alignment(combined_results, ref_cmaps, mol_maps, outname, argstring, multimap_mols);
//    combined_results.clear();
//    cout << "Writing partials\n";
//    if (partial_alignment) {
//        write_xmap_alignment(combined_results_partial, ref_cmaps, mol_maps, partialname, argstring, multimap_mols);
//    }
//    } else {
//    write_fda_alignment(combined_results, ref_cmaps, mol_maps, outname, multimap_mols);
//    if (partial_alignment) {
//        write_fda_alignment(combined_results_partial, ref_cmaps, mol_maps, partialname, multimap_mols);
//    }

//    chrono::steady_clock::time_point outWallE = chrono::steady_clock::now();

    double readWall = chrono::duration_cast<chrono::milliseconds>(readWallE - readWallS).count()/1000.;
    double ppWall = chrono::duration_cast<chrono::milliseconds>(ppWallE - ppWallS).count()/1000.;
    double alnWall = chrono::duration_cast<chrono::milliseconds>(alnWallE - alnWallS).count()/1000.;
//    double outWall = chrono::duration_cast<chrono::milliseconds>(outWallE - outWallS).count()/1000.;

    clock_t end = clock();
    chrono::steady_clock::time_point endWall = chrono::steady_clock::now();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double elapsedWall = chrono::duration_cast<chrono::milliseconds>(endWall - beginWall).count();
    printf("Total multithreaded CPU time for run: %.3f seconds\n",elapsed_secs);
    printf("Total CPU wall time for run: %.3f seconds\n",elapsedWall);
    cout << "Estimated wall time breakdown (seconds): \n";
    printf("Reading data: %.3f\n",readWall);
    printf("Preprocessing data: %.3f\n",ppWall);
    printf("Generating seeds and alignments: %.3f\n",alnWall);
//    printf("Writing alignments: %.3f\n", outWall);
    cout << endl;

    return 0;
}