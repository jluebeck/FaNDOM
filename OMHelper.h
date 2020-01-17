using namespace std;

#ifndef FANDOM_OMHELPER_H
#define FANDOM_OMHELPER_H

vector<pair<int,int>> mergeIntervals(vector<pair<int,int>> ivect);

map<int,vector<seedData>> mol_seeds_to_aln_regions(vector<seedData> &mol_seeds,
        map<int,vector<double>> &ref_cmaps, vector<double> &mol_posns, double padding);

void filter_mols(map<int,vector<double>> &mol_map, int min_map_lab, int min_map_len);

//Makes reverse cmap_map entries and adds them into the map
map<int,vector<double>> make_reverse_cmap(map<int,vector<double>> &cmap_map);

#endif //FANDOM_OMHELPER_H



