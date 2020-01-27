#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <ctime>
#include <set>
#include <future>
#include <chrono>

#include "threadsafe_queue.h"
#include "OMIO.h"
#include "OMHelper.h"
#include "fandomAligner.h"

using namespace std;

const int lookback = 5;
const int min_map_lab = 10;
const int min_map_len = 25000;
//const double fact_lookup [11] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
int n_threads = 6;
int min_aln_len = 6;
bool subsect = false;
float aln_padding = 5000;
//array<float,16> bins;
map<int,vector<pair<int,int>>> bed_data;

//------------------------------------------------------------------------------------------
//Do bookkeeping and manage work for each thread
map<int, vector<Alignment>> run_aln(map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &mols,
        map<int, vector<seedData>> &mol_seed_data, threadsafe_queue<int> &mol_ids_queue) {

    map<int, vector<Alignment>> results;
    while(true) {
        int mol_id = mol_ids_queue.pop();
//        if (mol_ids_queue.empty()) {
        if (!mol_id) {
            break;
        } else {
//            const int mol_id = mol_ids_queue.pop();
            //DEBUG PRINT vvvv
            cout << mol_id << "\n";
            int qsize = mol_ids_queue.size();
            if (qsize % 10000 == 0) {
                cout << qsize << " molecules left in queue\n";
            }
            double best_score = 0;
            vector<Alignment> mol_alns;
            vector<double> mol_vect = mols[mol_id];
            vector<seedData> curr_mol_seed_data = mol_seed_data[mol_id];
            map<int, vector<seedData>> curr_mol_seeds = mol_seeds_to_aln_regions(curr_mol_seed_data, ref_cmaps, mol_vect, aln_padding);
            for (auto &x: curr_mol_seeds) {
                int ref_id = x.first;
                vector<seedData> seed_list = x.second;

                vector<double> ref_vect = ref_cmaps[ref_id];
                for (const auto &s: seed_list) {
                    int a = s.ref_aln_lb;
                    int b = s.ref_aln_rb;
//                    cout << mol_id << " " << ref_id << " " << a << "," << b << "\n";
                    vector<vector<double>> S((b - a), vector<double>(mol_vect.size() - 1, 0));
                    vector<vector<pair<int, int>>> previous((b - a),
                                                            vector<pair<int, int>>(mol_vect.size() - 1, {-1, -1}));
                    dp_aln(S, previous, mol_vect, ref_vect, a, b, lookback);
                    Alignment curr_aln = dp_backtracking(S, previous, a, ref_id, mol_id);
                    curr_aln.seed_num = s.seed_num;
//                    cout << "f\n";

                    //check alignment length to see if usable
//                    cout << curr_aln.alignment.size() << " asize " << curr_aln.ref_id << " " << curr_aln.mol_id << " " << a << " " << b << " " << mol_vect.back() << "\n";
                    if (curr_aln.alignment.empty()) {
                        cout << "0 len alignment for " << curr_aln.ref_id << " " << curr_aln.mol_id << " " << a << " " << b << " " << mol_vect.back() << "\n";
                        continue;
                    }
                    double curr_score = get<2>(curr_aln.alignment.back());
                    if (curr_aln.alignment.size() < min_aln_len || curr_score/curr_aln.alignment.size() < 5000) { //TODO: ADD A BETTER SCORE THRESHOLD
                        continue;
                    }
//                    cout << "f1\n";

                    if (curr_score > best_score) {
                        best_score = curr_score;
                    }
                    mol_alns.push_back(curr_aln);
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
    }
//    prom_obj->set_value(results);
    return results;
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
    string ref_cmap_file, queryfile, bedfile, keyfile, sample_name, seed_file;

    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]).rfind("-t=", 0) == 0) {
            n_threads = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

        } else if (string(argv[i]).rfind("-r=", 0) == 0) {
            ref_cmap_file = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-q=", 0) == 0) {
            queryfile = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-b=", 0) == 0) {
            bedfile = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-k=", 0) == 0) {
            keyfile = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-sname=", 0) == 0) {
            sample_name = string(argv[i]).substr(string(argv[i]).find('=') + 1);

        } else if (string(argv[i]).rfind("-seeds=", 0) == 0) {
            seed_file = string(argv[i]).substr(string(argv[i]).find('=') + 1);
        }
    }
    return make_tuple(ref_cmap_file,queryfile,bedfile,keyfile,sample_name,seed_file);
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
    string ref_cmap_file, queryfile, bedfile, keyfile, sample_name, seed_file;
    tie(ref_cmap_file,queryfile,bedfile,keyfile, sample_name, seed_file) = parse_args(argc,argv);
    if (!bedfile.empty()) {
        subsect = true;
    }

    if (ref_cmap_file.empty()) {
        cout << "Reference genome .cmap file [-r=] unspecified\n";
        exit(0);
    } else if (queryfile.empty()) {
        cout << "Molecule .bnx file [-b=] unspecified\n";
        exit(0);
    } else if (sample_name.empty()) {
        cout << "Sample name [-s=] unspecified\n";
    } else if (seed_file.empty()) {
        cout << "Seed file must be specified for now\n";
        exit(0);
    } else if (bedfile.empty()) {
        cout << "Bed file not supplied, not subsecting alignments\n";
    }

    for (const string &fname: {ref_cmap_file,queryfile,seed_file}) {
        if (!file_exists(fname)) {
            cout << fname << " does not exist\n";
            exit(1);
        }
    }

    if (keyfile.empty()) {
        string rcf_base = ref_cmap_file.substr(ref_cmap_file.find('.'),ref_cmap_file.length());
        keyfile = rcf_base + "_key.txt";
    }

    cout << "reference genome cmap file: " << ref_cmap_file << "\n";
    cout << "query file: " << queryfile << "\n";
    cout << "\nNumber of threads for alignment: " << n_threads << "\n";
    if (bedfile.empty()) {
        cout << "no bed file supplied for subsecting\n";
    } else {
        cout << "bed file: " << bedfile << "\n" << "key file: " << keyfile << "\n";
    }

    //------------------------
    //parse inputs
    //parse ref genome
    chrono::steady_clock::time_point readWallS = chrono::steady_clock::now();

    map<int,vector<double>> ref_genome_unrev = parse_cmap(ref_cmap_file);
    cout << "Finished reading the ref\n";
    map<int,vector<double>> ref_cmaps = make_reverse_cmap(ref_genome_unrev);
//    for (auto &x: ref_cmaps) {
//        cout << x.first << " " << x.second.size() << "\n";
//    }

    //parse key and bed
    map<string,pair<int,double>> key_data = parse_key(keyfile);
    bed_data = parse_bed(bedfile,key_data);

    //parse seeds
    map<int, vector<seedData>> mol_seed_data = parse_seeds(seed_file,ref_cmaps);
    cout << "Read seeds from " << mol_seed_data.size() << " molecules\n";

    //parse queryfile
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
    //Preprocessing for alignments

    //make the mol_id_queue
    chrono::steady_clock::time_point ppWallS = chrono::steady_clock::now();
    filter_mols(mol_maps,min_map_lab,min_map_len);
    cout << "Got " << mol_maps.size() << " molecules passing initial filters\n";

    threadsafe_queue<int> mol_id_queue;
    //int test_count = 0;
    for (const auto &i: mol_maps) {
        int mol_id = i.first;
        //CHECK IF IT'S IN THE SEEDS
        if (mol_seed_data.find(mol_id) != mol_seed_data.end()) {
            mol_id_queue.push(mol_id);
//            break;
        }
        //DEBUG
//        if (test_count > 500) {
//            break;
//        }
//        test_count++;
    }

    //TODO: REMOVE ALL SEEDS OUTSIDE BED DATA
    chrono::steady_clock::time_point ppWallE = chrono::steady_clock::now();

    //------------------------------------------------------
    //launch alignments
    chrono::steady_clock::time_point alnWallS = chrono::steady_clock::now();
    cout << "Performing alignments\n";
    vector<future< map<int, vector<Alignment>>>> futs;
    vector<promise< map<int, vector<Alignment>>>> promises(n_threads);
    for (int i = 0; i < n_threads; i++) {
        futs.push_back(async(launch::async, run_aln, ref(ref_cmaps), ref(mol_maps), ref(mol_seed_data), ref(mol_id_queue)));
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
    cout << "Finished molecule alignment. \n";
    chrono::steady_clock::time_point alnWallE = chrono::steady_clock::now();

    //------------------------------------------------------
    //write output
    //TODO: write this on the fly
    chrono::steady_clock::time_point outWallS = chrono::steady_clock::now();
    string outname = sample_name + "_aln.txt";
    cout << "Writing alignments\n";
    write_alignment(combined_results,ref_cmaps,mol_maps,outname);
    chrono::steady_clock::time_point outWallE = chrono::steady_clock::now();

    //DEBUG PRINT:
//    for (auto i: mol_aln_windows) {
//        int mol_id = i.first;
//        cout << mol_id << "\n";
//        map<int,vector<pair<int,int>>> aln_windows = i.second;
//        for (auto k: aln_windows) {
//            cout << k.first << "\n";
//            for (auto l: k.second) {
//                cout << l.first << "-" << l.second << " ";
//            }
//            cout << "\n";
//        }
//        cout << "\n";
//    }

    double readWall = chrono::duration_cast<chrono::milliseconds>(readWallE - readWallS).count()/1000.;
    double ppWall = chrono::duration_cast<chrono::milliseconds>(ppWallE - ppWallS).count()/1000.;
    double alnWall = chrono::duration_cast<chrono::milliseconds>(alnWallE - alnWallS).count()/1000.;
    double outWall = chrono::duration_cast<chrono::milliseconds>(outWallE - outWallS).count()/1000.;

    clock_t end = clock();
    chrono::steady_clock::time_point endWall = chrono::steady_clock::now();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double elapsedWall = chrono::duration_cast<chrono::seconds>(endWall - beginWall).count();
    printf("Total multithreaded CPU time for run: %.3f seconds\n",elapsed_secs);
    printf("Total CPU wall time for run: %.3f seconds\n",elapsedWall);
    cout << "Estimated wall time breakdown (seconds): \n";
    printf("Reading data: %.3f\n",readWall);
    printf("Preprocessing data: %.3f\n",ppWall);
    printf("Generating alignments (includes seed boundary identification): %.3f\n",alnWall);
    printf("Writing alignments: %.3f\n", outWall);
    cout << endl;

    return 0;
}