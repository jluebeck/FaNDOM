#include "fandomAligner.h"

using namespace std;
float dist_scale_par=1.28;
int penalty_par=6000;

//fast pow function written by Martin Ankerl
//https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPrecisePow(double a, double b) {
    // calculate approximation with fraction of the exponent
    int e = (int) b;
    union {
        double d;
        int x[2];
    } u = { a };
    u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;

    // exponentiation by squaring with the exponent's integer part
    // double r = u.d makes everything much slower, not sure why
    double r = 1.0;
    while (e) {
        if (e & 1) {
            r *= a;
        }
        a *= a;
        e >>= 1;
    }

    return r * u.d;
}

//Alignment
double basic_score(double d_r, double d_m, int u_m, int u_r) {
//    double delta = pow(abs(d_r - d_m), dist_scale_par);
    double delta = fastPrecisePow(abs(d_r - d_m), dist_scale_par);
    double fl = penalty_par*(u_r + u_m);
    return 10000 - (fl + delta);
}

////initialize vector
//vector<vector<double>> S((b-a),vector<double> (mol_vect.size()-1,0));
//vector<vector<pair<int,int>>> previous((b-a),vector<pair<int,int>> (mol_vect.size()-1,{-1,-1}));
void dp_aln(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, vector<double> &mol_vect,
            vector<double> &ref_vect, int a, int b, int lookback) {

for (int j_ind = 0; j_ind < (b-a); j_ind++) {
    int i_start = max(0,j_ind-lookback);
        for (int q_ind = 0; q_ind < mol_vect.size()-1; q_ind++) { //Assumes last position is size
            int p_start = max(0,q_ind-lookback);
            for (int i_ind = i_start; i_ind < j_ind; i_ind ++) {
                double d_r = ref_vect[j_ind+a] - ref_vect[i_ind+a];
                int u_r = (j_ind - i_ind) - 1;
                for (int p_ind = p_start; p_ind < q_ind; p_ind++) {
                    double d_m = mol_vect[q_ind] - mol_vect[p_ind];
                    int u_m = (q_ind - p_ind) - 1;
                    double curr_score = S[i_ind][p_ind] + basic_score(d_m,d_r,u_m,u_r);
                    if (curr_score > S[j_ind][q_ind]) {
                        S[j_ind][q_ind] = curr_score;
                        previous[j_ind][q_ind] = {i_ind,p_ind};
                    }
                }
            }
        }
    }
}

pair<int,int> get_max_pair(const vector<vector<double>> &S, int ref_width, int mol_len) {
    double max_score = -100000;
    pair<int,int> max_pair = {-1,-1};
    for (int j_ind = 0; j_ind < ref_width; j_ind++) {
        for (int q_ind = 0; q_ind < mol_len; q_ind++) {
            if (S[j_ind][q_ind] > max_score) {
                max_score = S[j_ind][q_ind];
                max_pair = {j_ind,q_ind};
            }
        }
    }

    return max_pair;
}

Alignment dp_backtracking(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, pair<int,int> max_pair,
                          int a, int ref_id, int mol_id) {
    pair<int,int> null_tup = {-1,-1};
    //    pair<int,int> max_tup = get_max_pair(S);

    //backtracking
    vector<tuple<int,int,long>> aln;
    pair<int,int> curr_tup = max_pair;
    while (curr_tup != null_tup) {
        int j_ind = curr_tup.first;
        int q_ind = curr_tup.second;
        double score = S[j_ind][q_ind];
        aln.emplace_back(j_ind+a,q_ind,lround(score));
        curr_tup = previous[j_ind][q_ind];
    }
    reverse(aln.begin(),aln.end());

    Alignment align_struct(aln,ref_id,mol_id);
    return align_struct;
}