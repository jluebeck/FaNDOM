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
int score_limit = 5000;
//Data filtering
int Partial_aligning = 0;
int n_threads = 1;
int min_aln_len = 6;
int partial_alignment = 0;
float aln_prop_thresh_to_remap = 0.8;
double aln_len_thresh_to_remap = 25000.;
//bool subsect = false;
map<int,vector<pair<int,int>>> bed_data;
ofstream myfile;

//------------------------------------------------------------------------------------------
//Run Fandom alignment on a batch of seeds
map<int, vector<Alignment>> run_aln(map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &mols,
        map<int, vector<seedData>> &mol_seed_data) {
    map<int, vector<Alignment>> results;
    for (auto &item: mol_seed_data) {
        int mol_id = item.first;
        vector<seedData> curr_mol_seed_data = item.second;
        double best_score = 0;
        vector<Alignment> mol_alns;
        vector<double> mol_vect = mols[mol_id];
        //TODO: Siavash, this following function in OMHelper creates an alignment window from the seeds it is given. This step can be optimized
        /// Siavash: How??
        // check lb and rb are not zero
        map<int, vector<seedData>> curr_mol_seeds = mol_seeds_to_aln_regions(curr_mol_seed_data, ref_cmaps, mol_vect, aln_padding);
        for (auto &x: curr_mol_seeds) {
            int ref_id = x.first;
            vector<seedData> seed_list = x.second;
            vector<double> ref_vect = ref_cmaps[ref_id];
            for (const auto &s: seed_list) {
                int a = s.ref_aln_lb;
                int b = s.ref_aln_rb;
                vector<vector<double>> S((b - a), vector<double>(mol_vect.size() - 1, 0));
                vector<vector<pair<int, int>>> previous((b - a),
                                                        vector<pair<int, int>>(mol_vect.size() - 1, {-1, -1}));
                dp_aln(S, previous, mol_vect, ref_vect, a, b, lookback);
                pair<int,int> max_pair =  get_max_pair(S);
                double curr_score = S[max_pair.first][max_pair.second];
                if (multimap_mols || curr_score > best_score || partial_alignment) {
                    Alignment curr_aln = dp_backtracking(S, previous, max_pair, a, ref_id, mol_id);
                    curr_aln.seed_num = s.seed_num;
                    if (partial_alignment) {
                        curr_aln.is_partial = 1;
                    }
                    if (curr_aln.alignment.size() < min_aln_len ||
                        curr_score / curr_aln.alignment.size() < score_limit) { //TODO: ADD A BETTER SCORE THRESHOLD
                        continue;
                    }
                    mol_alns.push_back(curr_aln);
                }
                if (curr_score > best_score) {
                    best_score = curr_score;
                }
            }
        }
        //check if secondary alignments and set the status
        int is_multimapped = (int) (mol_alns.size() > 1);
        for (auto &curr_aln: mol_alns) {
            curr_aln.is_multimapped = is_multimapped;
            if (get<2>(curr_aln.alignment.back()) < best_score) {
                curr_aln.is_secondary = 1;
            }
        }
        results[mol_id] = mol_alns;
    }
    return results;
}

//run the filter and call the alignment
map<int, vector<Alignment>> filt_and_aln(int thread_num, map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &mols,
                                        map<int, dis_to_index> &ref_DTI, map<int, int> &ref_lens,
                                        threadsafe_queue<int> &mol_id_queue) {
    int counter = 0;
    vector<query> qq;
    map<int, vector<seedData>> seed_batch;
    map<int, vector<Alignment>> curr;
    map<int, vector<Alignment>> result;
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
        if (counter > 5000) {
            seed_batch = OMFilter(qq, thread_num, ref_DTI, ref_lens, ref_cmaps, mols, partial_alignment);
            curr = run_aln(ref_cmaps, mols, seed_batch);
            //////Log
            for (auto key : curr){
                myfile<< key.first<<"\n";
            }
            //////
            result.insert(curr.begin(), curr.end());
            counter = 0;
            qq.clear();
        }
    }
    //cleanup
    if (! qq.empty()) {
        seed_batch = OMFilter(qq, thread_num, ref_DTI, ref_lens, ref_cmaps, mols, partial_alignment);
        curr = run_aln(ref_cmaps, mols, seed_batch);
        //////Log
        for (auto key : curr){
            cout<< key.first<<endl;
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
//handle argss
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
        } else if ((string(argv[i]).rfind("-PA=", 0) == 0)) {
            Partial_aligning = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
        } else if ((string(argv[i]).rfind("-outfmt=", 0) == 0)) {
            outfmt = string(argv[i]).substr(string(argv[i]).find('=') + 1);
        } else if ((string(argv[i]).rfind("-version", 0) == 0)) {
            cout << "FaNDOM version " << version << "\n";
            exit(0);
        } else if ((string(argv[i]).rfind("-multimap", 0) == 0)) {
            multimap_mols = true;
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

    myfile.open(log_dir);
    myfile<<"KIR"<<endl;
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

    chrono::steady_clock::time_point ppWallE = chrono::steady_clock::now();

    //------------------------------------------------------
    //launch alignments
    chrono::steady_clock::time_point alnWallS = chrono::steady_clock::now();
    cout << "Performing alignments\n";
    vector<future< map<int, vector<Alignment>>>> futs;
    vector<promise< map<int, vector<Alignment>>>> promises(n_threads);
   for (int i = 0; i < n_threads; i++) {
       futs.push_back(async(launch:: async, filt_and_aln, i, ref(ref_cmaps), ref(mol_maps), ref(ref_DTI),
               ref(ref_num_to_length), ref(mol_id_queue)));
   }
    //-------------------------------------------------------
    //gather results
    vector<Alignment> combined_results;
   for (auto &f: futs) {
       map<int, vector<Alignment>> curr_result = f.get();
       for (const auto &x: curr_result) {
           combined_results.reserve(combined_results.size() + distance(x.second.begin(),x.second.end()));
           combined_results.insert(combined_results.end(),x.second.begin(),x.second.end());
       }
   }
    chrono::steady_clock::time_point alnWallE = chrono::steady_clock::now();
    cout << "Finished molecule alignment. \n" << combined_results.size() << " total alignments\n";
    vector<Alignment> combined_results_partial;
    //------------------------------------------------------
    if(Partial_aligning==1){
        partial_alignment = 1;
        unordered_set<int> mols_to_remap = get_remap_mol_ids(combined_results, mol_maps, aln_prop_thresh_to_remap,
                                                             aln_len_thresh_to_remap);
        cout << mols_to_remap.size() << " molecules will undergo partial-seeding.\n";
        myfile << mols_to_remap.size() << " molecules will undergo partial-seeding.\n";
        ////////////////////
        score_limit =1000;
        re_scale_par = 1.4;
        penalty_par = 7500;
        //Here again make thread and run partial alignments for remaining molecules
        if (mols_to_remap.size() > 0) {
            cout << "Performing partial alignments\n";
            myfile<< "Performing partial alignments\n";
            threadsafe_queue<int> mol_id_queue_partial;
            for (auto &i: mols_to_remap) {
                mol_id_queue_partial.push(i);
            }
            vector<future< map<int, vector<Alignment>>>> futs_partial;
            vector<promise< map<int, vector<Alignment>>>> promises_partial(n_threads);
            for (int i = 0; i < n_threads; i++) {
                futs_partial.push_back(async(launch::async, filt_and_aln, i, ref(ref_cmaps), ref(mol_maps), ref(ref_DTI),
                                             ref(ref_num_to_length), ref(mol_id_queue_partial)));
            }
            for (auto &f: futs_partial) {
                map<int, vector<Alignment>> curr_result = f.get();
                for (const auto &x: curr_result) {
                    combined_results_partial.reserve(combined_results_partial.size() + distance(x.second.begin(), x.second.end()));
                    combined_results_partial.insert(combined_results_partial.end(), x.second.begin(), x.second.end());
                }
            }
        }
    }
    ////////////////////
    //write output
    //TODO: write this on the fly
    chrono::steady_clock::time_point outWallS = chrono::steady_clock::now();
    string outname = sample_name;
    cout << "Writing alignments\n";
    if (outfmt == "xmap") {
        string fulname = "";
        fulname= outname+ ".xmap";
        string argstring;
        for (int i = 0; i < argc; ++i) {
            argstring += " ";
            argstring += argv[i];
        }
        write_xmap_alignment(combined_results, ref_cmaps, mol_maps, fulname, argstring, multimap_mols);
        if (Partial_aligning ==1 ) {
            string partial_name = outname + "_partial.xmap";
            write_xmap_alignment(combined_results_partial, ref_cmaps, mol_maps, partial_name, argstring, multimap_mols);
        }
    } else {
        string fulname = "";
        fulname= outname+ ".fda";
        write_fda_alignment(combined_results, ref_cmaps, mol_maps, fulname, multimap_mols);
        if (Partial_aligning ==1 ) {
            string partial_name = outname + "_partial.fda";
            write_fda_alignment(combined_results_partial, ref_cmaps, mol_maps, partial_name, multimap_mols);
        }
    }
    chrono::steady_clock::time_point outWallE = chrono::steady_clock::now();

    double readWall = chrono::duration_cast<chrono::milliseconds>(readWallE - readWallS).count()/1000.;
    double ppWall = chrono::duration_cast<chrono::milliseconds>(ppWallE - ppWallS).count()/1000.;
    double alnWall = chrono::duration_cast<chrono::milliseconds>(alnWallE - alnWallS).count()/1000.;
    double outWall = chrono::duration_cast<chrono::milliseconds>(outWallE - outWallS).count()/1000.;

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
    printf("Writing alignments: %.3f\n", outWall);
    cout << endl;

    return 0;
}