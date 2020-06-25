//AUTHOR: Siavash Raisi, sraeisid@ucsd.edu

#ifndef FANDOM_OMFILTER_H
#define FANDOM_OMFILTER_H

#include <vector>
#include <map>
#include <queue>
#include <cmath>
//#include <chrono>
#include <algorithm>

#include "SeedGraph.h"
#include "OMIO.h"

using namespace std;

extern int w;
extern int tolerance;
extern int ranked;
extern int ranked_partial;
extern int threshold;
extern int band_width;

///////All Definiation of typedef//////
typedef pair<int, pair<int, int>> score_ref_pos;
typedef vector<pair<int, int>> vec_pair_int;
typedef map<double, vector<pair<int, int>>> dis_to_index;
typedef vector<vector<pair<int, int>>> v_v_p;
//////////////////////////////////////////////

////////Structs
// This struct is for each possible candidate for doing alignments contain score of locations, refrence contig, direction, position and chain
struct answer {
    long score;
    int ref_contig;
    int pos;
    char dir;
//    vector<pair<int, int>> pairs;

    bool operator<(const answer &rhs) const {
        return score < rhs.score;
    }
};

//This struct is for saving possible break point which has reference contig, direction, position in terms of basepair, label position and score
struct breakpoint {
    int ref_contig;//chromosome
    char dir;
    double pos;
    int label_pos;
    float score;
    float long_score;
    int f1;//refrence
    int e1;//refrence
    int query_start;
    int query_end;
    bool operator<(const breakpoint &rhs) const {
        return long_score > rhs.long_score;
    }
};
bool compareBP(breakpoint a , breakpoint b );
//This struct is  for each query
struct query {
    int number; // query id
    dis_to_index distance; // genomic distance between pairs
    int length;// query length
    int start_pont; // In the big table, whats is the start point of this query
    v_v_p seeds_straight;// each refrence label has a label of pain which is seeds in straight direction. Vector of position in ref to vector of seeds
    v_v_p seeds_reverse;// each refrence label has a label of pain which is seeds in reverse direction. Vector of position in ref to vector of seeds
    priority_queue<answer> ans;// heap for saving top best candidate locations for doing alignment
    int find_loc;// boolean for show if we find great location to stop searching for this query
//    map<int, vector<breakpoint>> bp_from_start;// for each label possible break point from start to this label. Map label to vector of breakpint
//    map<int, vector<breakpoint>> bp_to_end;// for each label possible break point from label to end of query. Map label to vector of breakpint
    priority_queue<breakpoint> bp;
    query(int n, dis_to_index d, int l , int find){
        number = n;
        distance = d;
        length = l;
        find_loc = find;
    }
};

///////////////OMFilter Functions
map<double, vec_pair_int> calculate_genomic_distance(vector<double> dist);

map<int, dis_to_index> genomic_distance(map<int, vector<double>> &gen);

int **allocateTwoDimenArrayOnHeapUsingMalloc(int row, int col);

void destroyTwoDimenArrayOnHeapUsingFree(int **ptr, int row, int col);

int **merge_list(dis_to_index LM, dis_to_index LN);

//void calculate_seeds_score_in_band(vector<query> &qq, int ref_contig);
void calculate_seeds_score_in_band(vector<query> &qq, int ref_contig, vector<double> &ref_dist,
        map<int, vector<double>> &query_genome);

//void calculate_seeds_score_in_band_SV(vector<query> &qq, int ref_contig) {
void calculate_seeds_score_in_band_SV(vector<query> &qq, int ref_contig, vector<double> &ref_dist,
        map<int, vector<double>> &query_genome);

map<int, vector<seedData>> OMFilter(vector<query> qq, int number, map<int, dis_to_index> &ref_list,
map<int, int> &ref_lens, map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &query_cmaps,
int SV_detection);

#endif //FANDOM_OMFILTER_H
