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
#include <unistd.h>
#include <cmath>

#include "OMHelper.h"

using namespace std;

//inline bool file_exists (const string &fname) {
//    ifstream f(fname.c_str());
//    return f.good();
//}


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

vector<vector<float>> parse_bed(const string &fname);
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

void write_xmap_header(ofstream &outfile, const string &argstring);

void write_fda_header(ofstream &outfile);

void write_xmap_alignment(ofstream &outfile, vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref,
        map<int,vector<double>> &mol_map, bool is_mm, size_t a_pair_ind);

void write_fda_alignment(ofstream &outfile, vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref,
        map<int,vector<double>> &mol_map, bool is_mm, size_t a_pair_ind);

#endif //FANDOM_OMHELPER_H
