//AUTHOR: Jens Luebeck, jluebeck@ucsd.edu

#ifndef FANDOM_OMHELPER_H
#define FANDOM_OMHELPER_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_set>

using namespace std;

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

struct Alignment {
    vector<tuple<int,int,double>> alignment;
    int ref_id;
    int mol_id;
    int is_multimapped = 0;
    int is_secondary = 0;
    int seed_num = 0;

    Alignment(vector<tuple<int,int,double>> &new_alignment, int R_ID, int M_ID) {
        alignment = new_alignment; ref_id = R_ID; mol_id = M_ID;
    }
};

vector<pair<int,int>> mergeIntervals(vector<pair<int,int>> ivect);

map<int,vector<seedData>> mol_seeds_to_aln_regions(vector<seedData> &mol_seeds,
        map<int,vector<double>> &ref_cmaps, vector<double> &mol_posns, double padding);

void filter_mols(map<int,vector<double>> &mol_map, int min_map_lab, int min_map_len);

//get molecule ids needing a partial alignment seeding
unordered_set<int> get_remap_mol_ids(const vector<Alignment> &aln_list, map<int, vector<double>> &mol_maps,
                              float remap_prop_cut, double remap_len_cut);

//Makes reverse cmap_map entries and adds them into the map
map<int,vector<double>> make_reverse_cmap(map<int,vector<double>> &cmap_map);

map<int, int> calculate_length(const map<int, vector<double>> &contigs);

#endif //FANDOM_OMHELPER_H



