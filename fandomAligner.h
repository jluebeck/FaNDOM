//AUTHOR: Jens Luebeck, jluebeck@ucsd.edu

#ifndef FANDOM_FANDOMALIGNER_H
#define FANDOM_FANDOMALIGNER_H

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <map>
#include <fstream>
//#include "OMIO.h"
#include "OMHelper.h"

using namespace std;

//
double basic_score(double d_r, double d_m, int u_m, int u_r);

//
void dp_aln(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, vector<double> &mol_vect,
            vector<double> &ref_vect, int a, int b, int lookback);

//
Alignment dp_backtracking(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, int a, int ref_id, int mol_id);

#endif //FANDOM_FANDOMALIGNER_H
