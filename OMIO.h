//AUTHOR: Jens Luebeck, jluebeck@ucsd.edu

#ifndef FANDOM_OMIO_H
#define FANDOM_OMIO_H

#include <iostream>
#include <cstdio>
#include <fstream>

#include <vector>
#include <map>
#include <ctime>
#include <set>
#include <sstream>

using namespace std;

//inline bool file_exists (const string &fname) {
//    ifstream f(fname.c_str());
//    return f.good();
//}

struct seedData {
    int ref_id;
    int mol_id;
    int ref_0_lab;
    int seed_num;
    int ref_aln_lb = 0;
    int ref_aln_rb = 0;
    seedData(int R_ID,int M_ID, int R_0_L, int SN) {
        ref_id = R_ID, mol_id = M_ID, ref_0_lab = R_0_L, seed_num = SN;
    }
};


inline bool file_exists (const string &fname) {
    ifstream ifile(fname.c_str());
    cout << fname << " " << ifile.good() << "\n";
    return ifile.good();
}

/**
 * Opens and parses a CMAP text file, stores cmap as pair which is given as input
 * @param fname
 * @param cmap_pairs
 */
map<int,vector<double>> parse_cmap(const string &fname);

///**
// * Prints a cmap map by map-id to cout
// */
void cmap_map_to_string(map<int,vector<double>> &cmap_map);

////Parse the score thresholds
//map<int,float> parse_score_thresholds(string &score_thresh_file);

//Parse list of contigs to aln with
set<int> parse_contig_list(string &contig_list_file);

//map<int,vector<pair<int,int>>> parse_bed(const string &fname, map<string,pair<int,float>> &key_data);
map<int,vector<pair<int,int>>> parse_bed(const string &fname, map<string,pair<int,double>> &key_data);

map<string,pair<int,double>> parse_key(const string &fname);

map<int,vector<double>> parse_bnx(const string &fname);

map<int, vector<seedData>> parse_seeds(const string &fname, map<int, vector<double>> &ref_cmaps);

#endif //FANDOM_OMHELPER_H
