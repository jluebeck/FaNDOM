#include "OMIO.h"

using namespace std;

#ifndef FANDOM_FANDOMALIGNER_H
#define FANDOM_FANDOMALIGNER_H
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

//
double basic_score(double d_r, double d_m, int u_m, int u_r);

//
void dp_aln(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, vector<double> &mol_vect,
            vector<double> &ref_vect, int a, int b, int lookback);

//
Alignment dp_backtracking(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, int a, int ref_id, int mol_id);

//
void write_alignment(vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref, map<int,vector<double>> &mol_map,
                     string &outname);

#endif //FANDOM_FANDOMALIGNER_H
