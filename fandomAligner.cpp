#include "fandomAligner.h"

using namespace std;

//Alignment
double basic_score(double d_r, double d_m, int u_m, int u_r) {
    double delta = pow(abs(d_r - d_m),1.1);
    double fl = 5000*(u_r + u_m);
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

Alignment dp_backtracking(vector<vector<double>> &S, vector<vector<pair<int,int>>> &previous, int a, int ref_id, int mol_id) {
    //get max element
    pair<int,int> null_tup = {-1,-1};
    double max_score = 0;
    pair<int,int> max_tup = {-1,-1};
    for (int j_ind = 0; j_ind < S.size(); j_ind++) {
        for (int q_ind = 0; q_ind < S[0].size(); q_ind++) {
            if (S[j_ind][q_ind] > max_score) {
                max_score = S[j_ind][q_ind];
                max_tup = {j_ind,q_ind};
            }
        }
    }

    //backtracking
    vector<tuple<int,int,double>> aln;
    pair<int,int> curr_tup = max_tup;
    while (curr_tup != null_tup) {
        int j_ind = curr_tup.first;
        int q_ind = curr_tup.second;
        double score = S[j_ind][q_ind];
        aln.emplace_back(j_ind+a,q_ind,score);
        curr_tup = previous[j_ind][q_ind];
    }
    reverse(aln.begin(),aln.end());

    Alignment align_struct(aln,ref_id,mol_id);
    return align_struct;
}

//TODO: REFACTOR
//TODO: ADD XMAP WRITING
void write_alignment(vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref, map<int,vector<double>> &mol_map,
                     string &outname) {

    ofstream outfile;
    outfile.open(outname);

    outfile << "#0\tref_id\tmol_id\taln_direction\tref_start_pos\tref_end_pos\tmol_start_pos\tmol_end_pos\tmol_length\n";
    outfile << "#1\ttotal_score\tmean_score\tis_multimapped\tis_secondary\taln_seed_num\n";
    outfile << "#2\talignment [aln_index]:(ref_pos, mol_pos, mol_lab, score_delta)\n";
    outfile << "#3\tcigar [aln_index]:(delta_ref, delta_mol, mol_label_diff, delta_difference)\n";

    size_t a_pair_ind = 0;
    for (const auto &aln_struct: aln_list) {
        int raw_ref_id = aln_struct.ref_id;
        int ref_id = abs(raw_ref_id);
        int mol_id = aln_struct.mol_id;
        vector<double> ref_posns = cmaps_ref[ref_id];
        vector<double> mol_posns = mol_map[mol_id];
        vector<tuple<int,int,double>> alignment = aln_struct.alignment;

        outfile << "0\t" << ref_id << "\t" << mol_id << "\t";
        //compute ref start and end posns
//        float ref_start_pos = cmaps_ref[raw_ref_id][get<0>(alignment[0])];
//        float ref_end_pos = cmaps_ref[raw_ref_id][get<0>(alignment.back())];
        float ref_start_pos;
        float ref_end_pos;
        double score = get<2>(alignment.back());
        double mean_score = score/(alignment.size()-1);
        char direction;
        if (raw_ref_id < 0) {
            direction = '-';
//            //reverse the alignment
//            float temp_start = ref_start_pos;
//            ref_start_pos = cmaps_ref[raw_ref_id].back() - ref_end_pos;
//            ref_end_pos = cmaps_ref[raw_ref_id].back() - temp_start;
            int ref_start_label = (ref_posns.size() - get<0>(alignment[0])) - 2;
            int ref_end_label = (ref_posns.size() - get<0>(alignment.back())) - 2;
            ref_start_pos = ref_posns[ref_start_label];
            ref_end_pos = ref_posns[ref_end_label];

        } else {
            direction = '+';
            ref_start_pos = ref_posns[get<0>(alignment[0])];
            ref_end_pos = ref_posns[get<0>(alignment.back())];

        }
        float mol_start_pos = mol_posns[get<1>(alignment[0])];
        float mol_end_pos = mol_posns[get<1>(alignment.back())];


        //WRITE 0
        //continued from above
        outfile << direction << "\t" << (int)ref_start_pos << "\t" << (int)ref_end_pos << "\t" << mol_start_pos << "\t" << mol_end_pos << "\t" << mol_map[mol_id].back() << "\n";

        //WRITE 1
        //score, multimapped (0/1), is_secondary (0/1)
        outfile << "1\t" << score << "\t" << mean_score << "\t" << aln_struct.is_multimapped << "\t" << aln_struct.is_secondary << "\t" << aln_struct.seed_num << "\n";

        string outstring2 = "2";
        string outstring3 = "3";
        float prev_score = 0.0;
        float prev_ref_pos = ref_start_pos;
        float prev_mol_pos = mol_start_pos;
        int prev_mol_lab = get<1>(alignment[0]);
        int aln_pos = 0;
        for (const auto &aln_tup: alignment) {
            int ref_lab = get<0>(aln_tup);
            int mol_lab = get<1>(aln_tup);
            float curr_score = get<2>(aln_tup);
            float score_delta = curr_score - prev_score;
            if (raw_ref_id < 0) {
                //translate the reverse label number
                int n_labs = cmaps_ref[raw_ref_id].size();
                ref_lab = (n_labs - ref_lab) - 2;
            }
            float rlp = ref_posns[ref_lab];
            float mlp = mol_posns[mol_lab];
            float abs_ref_diff = fabs(prev_ref_pos - rlp);
            float mol_diff = mlp-prev_mol_pos;


            //outfile << "#2\talignment (ref_pos, mol_pos, mol_lab, score_delta)\n";
            //outfile << "#3\tcigar (delta_ref, delta_mol, difference)\n";
            //WRITE 2
            outstring2+=("\t" + to_string(aln_pos) + ":(" + to_string(int(rlp)) + ", " + to_string(int(mlp)) + ", " + to_string(mol_lab+1) + ", " + to_string(score_delta) + ")");

            //WRITE 3

            outstring3+=("\t" + to_string(aln_pos) + ":(" + to_string(int(abs_ref_diff)) + ", " + to_string(int(mol_diff)) + ", " + to_string(mol_lab - prev_mol_lab - 1) + ", " + to_string(mol_diff - abs_ref_diff) + ")");

            prev_score = curr_score;
            prev_ref_pos = rlp;
            prev_mol_pos = mlp;
            prev_mol_lab = mol_lab;
            aln_pos++;
        }
        outfile << outstring2 << "\n" << outstring3 << "\n";
        a_pair_ind++;
    }
    outfile << flush;
    outfile.close();
}