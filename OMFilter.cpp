#include "OMFilter.h"

using namespace std;

int w= 3;
int tolerance = 350;
int ranked = 300;
//int a[57000][9000] = {0};   ///Unused?
int threshold = 3;
int band_width = 6000;



//////////////// Function for returning genomic distances of a list/////////
map<double, vec_pair_int> calculate_genomic_distance(vector<double> dist) {
    map<double, vec_pair_int> ans;
    for (int it = 0; it != dist.size(); it++) {
        for (int j = it + 1; j < min(it + w, int(dist.size())); j++) {
            if (ans.find(double(roundf((dist[j] - dist[it]) * 1) / 1)) == ans.end()) {
                vec_pair_int ind;
                ind.push_back(make_pair(it + 1, j + 1));
                ans.insert({double(roundf((dist[j] - dist[it]) * 1) / 1), ind});
            } else {
                ans[float(roundf((dist[j] - dist[it]) * 1) / 1)].push_back(
                        make_pair(it + 1, j + 1));
            }
        }
    }
    return ans;
}

/////////////////Calculate genomic distance of map of list (the same as calculate_genomic_distance function just do it for all key in a map)
map<int, dis_to_index> genomic_distance(map<int, vector<double>> &gen) {
    map<int, dis_to_index> genomic_dist;
    for (auto & it : gen) {
        genomic_dist.insert({it.first, calculate_genomic_distance(it.second)});
    }
    return genomic_dist;
}

////////Make a big 2d array////////////
int **allocateTwoDimenArrayOnHeapUsingMalloc(int row, int col) {
    int **ptr = (int **) malloc(sizeof(int *) * row);
    for (int i = 0; i < row; i++) {
        ptr[i] = (int *) malloc(sizeof(int) * col);
    }
    return ptr;
}

///////Delete 2d array/////////////
void destroyTwoDimenArrayOnHeapUsingFree(int **ptr, int row, int col) {
    for (int i = 0; i < row; i++) {
        free(ptr[i]);
    }
    free(ptr);
}

////////////This function merge list LM and LN in the paper
int **merge_list(dis_to_index LM, dis_to_index LN) {
    int **a = allocateTwoDimenArrayOnHeapUsingMalloc(55000, 8500);
    for (int i = 0; i < 55000; i++) {
        for (int j = 0; j < 8500; j++) {
            a[i][j] = 0;
        }
    }
    dis_to_index::iterator i = LN.begin();
    while (i != LN.end()) {
        int down = i->first - tolerance;//down tolerance
        int up = i->first + tolerance;//up tolerance
        dis_to_index::iterator down_index = LM.lower_bound(down);//doing binary search
        dis_to_index::iterator up_index = LM.upper_bound(up);//doing binary search
        set<int> q_p;// query point
        set<int> r_p;// reference point
        if (down_index != up_index) {
            for (dis_to_index::iterator it = down_index; it != up_index; it++) {
                for (int j = 0; j < it->second.size(); j++) {
                    for (int start_point = max(1, it->second[j].second - w + 1);
                         start_point <= it->second[j].first; start_point++) {
                        q_p.insert(start_point);//add to all window that contains this genomic distance in query
                    }
                }
            }
            for (int j = 0; j < i->second.size(); j++) {
                for (int s_p = max(1, i->second[j].second - w + 1); s_p <= i->second[j].first; s_p++) {
                    r_p.insert(s_p);// add to all window that containg this window in reference
                }
            }
            for (set<int>::iterator j = q_p.begin(); j != q_p.end(); j++) {
                for (set<int>::iterator k = r_p.begin(); k != r_p.end(); k++) {
                    a[*k][*j]++; /// increasing number in 2d array
                }
            }
        }
        i++;
    }
    return a;
}

////////////////////For finding complete alignment score/////////////
//void calculate_seeds_score_in_band(vector<query> &qq, int ref_contig) {
void calculate_seeds_score_in_band(vector<query> &qq, int ref_contig, vector<double> &ref_dist,
        map<int, vector<double>> &query_genome) {
    //    vector<double> ref_dist = ref_genome_raw[ref_contig];
    for (auto & q : qq) {
        if (q.find_loc == 0) {
            vector<double> q_dist = query_genome[q.number];
            double query_len = (q_dist.back() - q_dist[0]) * 1.414;
            int best_straight_score = 10000000;
            int prev_index = -1;
            answer best_straight_ans = answer{best_straight_score, ref_contig, -1, '+'};
            int best_rev_score = 10000000;
            int prev_index_rev = -1;
            int straight_counter = 0;
            int rev_counter = 0;
            answer best_rev_ans = answer{best_straight_score, ref_contig, -1, '+'};
            int last_index_str = 0;
            int last_index_rev = 0;
            for (int cc = 1; cc < q.seeds_straight.size(); cc++) {
                if (ref_dist[cc]-ref_dist[last_index_str]>3000) {
                    if (q.seeds_straight[cc].size() > 3) {
                        last_index_str = cc;
                        long straight_score = 10000000;
                        pair<long, vector<pair<int, int>>> score_path = solve_graph_straight(q.seeds_straight[cc], q_dist,
                                                                                             ref_dist, w);
                        straight_score = score_path.first;
                        if (abs(prev_index - cc) > 1) {
                            if (straight_counter > 0) {
                                double check_score = 1.0 - double(best_straight_score / query_len);
                                if (check_score > 0.65) {
                                    q.find_loc = 1;
                                    q.ans = priority_queue<answer>();
                                    q.ans.push(best_straight_ans);
                                    break;
                                } else if (q.ans.size() < ranked) {
                                    q.ans.push(best_straight_ans);
                                } else {
                                    if (q.ans.top().score > best_straight_score) {
                                        q.ans.pop();
                                        q.ans.push(best_straight_ans);
                                    }
                                }
                            }
                            straight_counter = 1;
                            best_straight_score = straight_score;
                            best_straight_ans = answer{straight_score, ref_contig, cc + 1, '+', score_path.second};
                            prev_index = cc;
                        } else {
                            if (straight_score < best_straight_score) {
                                best_straight_score = straight_score;
                                best_straight_ans = answer{straight_score, ref_contig, cc + 1, '+', score_path.second};
                            }
                            prev_index = cc;
                        }

                    }
                }
                if(ref_dist[cc]-ref_dist[last_index_rev]>3000){
                    if (q.seeds_reverse[cc].size() > 3) {
                        last_index_rev = cc;
                        long reverese_score = 10000000;
                        pair<long, vector<pair<int, int>>> score_path = solve_graph_reverse(q.seeds_reverse[cc], q_dist,
                                                                                            ref_dist, w);
                        reverese_score = score_path.first;
                        if (abs(prev_index_rev - cc) > 1) {
                            if (rev_counter > 0) {
                                double check_score = 1.0 - double(best_rev_score / query_len);
                                if (check_score > 0.65) {
                                    q.find_loc = 1;
                                    q.ans = priority_queue<answer>();
                                    q.ans.push(best_rev_ans);
                                    break;
                                } else if (q.ans.size() < ranked) {
                                    q.ans.push(best_rev_ans);
                                } else {
                                    if (q.ans.top().score > best_rev_score) {
                                        q.ans.pop();
                                        q.ans.push(best_rev_ans);
                                    }
                                }
                            }
                            rev_counter = 1;
                            best_rev_score = reverese_score;
                            best_rev_ans = answer{best_rev_score, ref_contig, cc + 1, '-', score_path.second};
                            prev_index_rev = cc;

                        } else {
                            if (reverese_score < best_rev_score) {
                                best_rev_score = reverese_score;
                                best_rev_ans = answer{best_rev_score, ref_contig, cc + 1, '-', score_path.second};
                            }
                            prev_index_rev = cc;
                        }
                    }
                }
            }
            if (q.ans.size() < ranked && best_straight_ans.pos != -1 && q.find_loc==0) {
                int ali = 0;
                q.ans.push(best_straight_ans);
            } else if (! q.ans.empty() ){
                if (q.ans.top().score > best_straight_score && best_straight_ans.pos != -1) {
                    q.ans.pop();
                    q.ans.push(best_straight_ans);
                }
            }
            if (q.ans.size() < ranked && best_rev_ans.pos != -1 &&q.find_loc==0) {
                q.ans.push(best_rev_ans);
            } else if (! q.ans.empty() ) {
                if (q.ans.top().score > best_rev_score && best_rev_ans.pos != -1) {
                    q.ans.pop();
                    q.ans.push(best_rev_ans);
                }
            }
        }
    }
}

////////////////////For finding partial alignment score/////////////
//void calculate_seeds_score_in_band_SV(vector<query> &qq, int ref_contig) {
void calculate_seeds_score_in_band_SV(vector<query> &qq, int ref_contig, vector<double> &ref_dist,
        map<int, vector<double>> &query_genome) {
//    vector<double> ref_dist = ref_genome_raw[ref_contig];
    for (auto & q : qq) {
        vector<double> q_dist = query_genome[q.number];
        double SV_minimum_length = max(40000.0, 0.15 * (q_dist[q_dist.size() - 1] - q_dist[0]));
        for (int cc = 0; cc < q.seeds_straight.size(); cc++) {
            if (q.seeds_straight[cc].size() > 1) {
                vector<long> score_from_first = solve_graph_straight_SV(q.seeds_straight[cc], q_dist, ref_dist, 1, w);
                vector<float> scores;
                for (int i = 0; i < q.seeds_straight[cc].size(); ++i) {
                    float partial_score = score_from_first[(i + 1) * w] /
                                          ((q_dist[q.seeds_straight[cc][i].first + w - 1 - 1] - q_dist[0]) * 1.414);
                    scores.push_back(partial_score);
                }
                vector<float> scores2;
                vector<long> score_to_end = solve_graph_straight_SV(q.seeds_straight[cc], q_dist, ref_dist, 0, w);
                for (int i = 0; i < q.seeds_straight[cc].size(); ++i) {
                    float partial_score = score_to_end[i * w + 1] /
                                          ((q_dist[q_dist.size() - 1] - q_dist[q.seeds_straight[cc][i].first - 1]) *
                                           1.414);
                    scores2.push_back(partial_score);
                }
                for (int i = 0; i < q.seeds_straight[cc].size(); ++i) {
                    if (scores[i] <= scores2[i]) {
                        int f = 0;
                        float partial_score = scores[i];
                        for (int j = 0; j < q.seeds_straight[cc].size(); ++j) {
                            if (scores[j] <= partial_score &&
                                q.seeds_straight[cc][i].first <= q.seeds_straight[cc][j].first && j != i) {
                                f = 1;
                                break;
                            }
                        }
                        if (partial_score < 0.75 &&
                            (q_dist[q.seeds_straight[cc][i].first + w - 1 - 1] - q_dist[0]) >= SV_minimum_length &&
                            f == 0) {
                            breakpoint new_bp = breakpoint{ref_contig, '+',
                                                           ref_dist[q.seeds_straight[cc][i].second + w - 1 - 1],
                                                           q.seeds_straight[cc][i].second + w - 1, partial_score};
                            int bp_index = q.seeds_straight[cc][i].first + w - 1;
                            int find = 0;
                            for (int fs_i = q.bp_from_start[bp_index].size() - 1; fs_i >= 0; fs_i--) {
                                breakpoint last_bp = q.bp_from_start[bp_index][fs_i];
                                if (last_bp.ref_contig == new_bp.ref_contig && last_bp.dir == new_bp.dir &&
                                    abs(last_bp.pos - new_bp.pos) < 20000) {
                                    find = 1;
                                    if (last_bp.score >= new_bp.score) {
                                        q.bp_from_start[bp_index][fs_i] = new_bp;
                                    }
                                }
                            }
                            if (find == 0) {
                                q.bp_from_start[bp_index].push_back(new_bp);
                            }
                            ////begu darim
                        }
                    }
                }

                for (int i = 0; i < q.seeds_straight[cc].size(); ++i) {
                    if (scores2[i] <= scores[i]) {
                        int f = 0;
                        float partial_score = scores2[i];
                        for (int j = 0; j < q.seeds_straight[cc].size(); ++j) {
                            if (scores2[j] <= partial_score &&
                                q.seeds_straight[cc][j].first < q.seeds_straight[cc][i].first && j != i) {
                                f = 1;
                                break;
                            }
                        }
                        if (partial_score < 0.75 &&
                            (q_dist[q_dist.size() - 1] - q_dist[q.seeds_straight[cc][i].first - 1]) >=
                            SV_minimum_length &&
                            f == 0) {
                            breakpoint new_bp = breakpoint{ref_contig, '+',
                                                           ref_dist[q.seeds_straight[cc][i].second - 1],
                                                           q.seeds_straight[cc][i].second, partial_score};
                            int bp_index = q.seeds_straight[cc][i].first;
                            int find = 0;
                            for (int fs_i = q.bp_to_end[bp_index].size() - 1; fs_i >= 0; fs_i--) {
                                breakpoint last_bp = q.bp_to_end[bp_index][fs_i];
                                if (last_bp.ref_contig == new_bp.ref_contig && last_bp.dir == new_bp.dir &&
                                    abs(last_bp.pos - new_bp.pos) < 20000) {
                                    find = 1;
                                    if (last_bp.score >= new_bp.score) {
                                        q.bp_to_end[bp_index][fs_i] = new_bp;
                                    }
                                }
                            }
                            if (find == 0) {
                                q.bp_to_end[bp_index].push_back(new_bp);
                            }
////begu darim

                        }
                    }
                }
            }
            if (q.seeds_reverse[cc].size() > 1) {
                vector<long> score_to_end = solve_graph_reverse_SV(q.seeds_reverse[cc], q_dist, ref_dist, 1, w);
                vector<float> scores;
                for (int i = 0; i < q.seeds_reverse[cc].size(); i++) {
                    float partial_score = score_to_end[(i + 1) * w] /
                                          ((q_dist[q_dist.size() - 1] -
                                            q_dist[q.seeds_reverse[cc][i].first - w + 1 - 1]) * 1.414);
                    scores.push_back(partial_score);
                }
                vector<long> score_from_first = solve_graph_reverse_SV(q.seeds_reverse[cc], q_dist, ref_dist, 0, w);
                vector<float> scores2;
                for (int i = 0; i < q.seeds_reverse[cc].size(); ++i) {
                    float partial_score = score_from_first[i * w + 1] /
                                          ((q_dist[q.seeds_reverse[cc][i].first - 1] - q_dist[0]) * 1.414);
                    scores2.push_back(partial_score);
                }
                for (int i = 0; i < q.seeds_reverse[cc].size(); i++) {
                    if (scores[i] <= scores2[i]) {
                        float partial_score = scores[i];
                        int f = 0;
                        for (int j = 0; j < q.seeds_reverse[cc].size(); j++) {
                            if (scores[j] <= partial_score &&
                                q.seeds_reverse[cc][j].first < q.seeds_reverse[cc][i].first &&
                                j != i) {
                                f = 1;
                                break;
                            }
                        }
                        if (partial_score < 0.75 &&
                            (q_dist[q_dist.size() - 1] - q_dist[q.seeds_reverse[cc][i].first - w + 1 - 1]) >=
                            SV_minimum_length && f == 0) {
                            breakpoint new_bp = breakpoint{ref_contig, '-',
                                                           ref_dist[q.seeds_reverse[cc][i].second + w - 1 - 1],
                                                           q.seeds_reverse[cc][i].second + w - 1, partial_score};
                            int bp_index = q.seeds_reverse[cc][i].first - w + 1;
                            int find = 0;
                            for (int fs_i = q.bp_to_end[bp_index].size() - 1; fs_i >= 0; fs_i--) {
                                breakpoint last_bp = q.bp_to_end[bp_index][fs_i];
                                if (last_bp.ref_contig == new_bp.ref_contig && last_bp.dir == new_bp.dir &&
                                    abs(last_bp.pos - new_bp.pos) < 20000) {
                                    find = 1;
                                    if (last_bp.score >= new_bp.score) {
                                        q.bp_to_end[bp_index][fs_i] = new_bp;
                                    }
                                }
                            }
                            if (find == 0) {
                                q.bp_to_end[bp_index].push_back(new_bp);
                            }
                            ///begu darim
                        }
                    }
                }

                for (int i = 0; i < q.seeds_reverse[cc].size(); ++i) {
                    if (scores2[i] <= scores[i]) {
                        float partial_score = scores2[i];
                        int f = 0;
                        for (int j = 0; j < q.seeds_reverse[cc].size(); ++j) {
                            if (scores2[j] <= partial_score &&
                                q.seeds_reverse[cc][i].first < q.seeds_reverse[cc][j].first && j != i) {
                                f = 1;
                                break;
                            }
                        }

                        if (partial_score < 0.75 &&
                            (q_dist[q.seeds_reverse[cc][i].first - 1] - q_dist[0]) >= SV_minimum_length && f == 0) {
                            breakpoint new_bp = breakpoint{ref_contig, '-', ref_dist[q.seeds_reverse[cc][i].second - 1],
                                                           q.seeds_reverse[cc][i].second, partial_score};
                            int bp_index = q.seeds_reverse[cc][i].first;
                            int find = 0;
                            for (int fs_i = q.bp_from_start[bp_index].size() - 1; fs_i >= 0; fs_i--) {
                                breakpoint last_bp = q.bp_from_start[bp_index][fs_i];
                                if (last_bp.ref_contig == new_bp.ref_contig && last_bp.dir == new_bp.dir &&
                                    abs(last_bp.pos - new_bp.pos) < 20000) {
                                    find = 1;
                                    if (last_bp.score >= new_bp.score) {
                                        q.bp_from_start[bp_index][fs_i] = new_bp;
                                    }
                                }
                            }
                            if (find == 0) {
                                q.bp_from_start[bp_index].push_back(new_bp);
                            }
                            ///begu darim
                        }
                    }
                }
            }
        }
    }
}

////////////////////
//formerly called do_aligning
//qq is a query list
//number is same as query_id
//ref_lens is ref_num_to_length

map<int, vector<seedData>> OMFilter(vector<query> qq, int number, map<int, dis_to_index> &ref_list, map<int, int> &ref_lens,
              map<int,vector<double>> &ref_cmaps, map<int,vector<double>> &query_cmaps, const int SV_detection) {

    vector<int> index_to_query = vector<int>(8501);
    map<int, vector<seedData>> molSeedMap;
    dis_to_index LM;
    int counter = 0;
    for (auto & q : qq) {
        for (auto & j : q.distance) {
            for (pair<int, int> &p : j.second) {
                p = make_pair(p.first + counter, p.second + counter);
            }
            if (LM.count(j.first) == 0) {
                LM[j.first] = j.second;
            } else {
                LM[j.first].insert(LM[j.first].end(), j.second.begin(), j.second.end());
            }
        }
        q.start_pont = counter;
        for (int y = counter; y < counter + q.length + w; y++) {
            index_to_query[y] = q.number;
        }
        counter = counter + w + q.length;
    }
    for (auto & ref_iterator : ref_list) {
        int ref_id = ref_iterator.first;
        for (auto & i : qq) {
            i.seeds_straight = v_v_p(ref_lens[ref_id] + 50);
            i.seeds_reverse = v_v_p(ref_lens[ref_id] + 50);
        }
        dis_to_index LN = ref_iterator.second;
//        chrono::steady_clock::time_point b1 = chrono::steady_clock::now();
        int **a = merge_list(LM, LN);
//        chrono::steady_clock::time_point b2 = chrono::steady_clock::now();
//        double b3 = chrono::duration_cast<chrono::milliseconds>(b2 - b1).count() / 1000.;
        // printf("making matrix: %.3f\n", b3);
        vector<double> ref_dist = ref_cmaps[ref_id];
        int l = ref_lens[ref_id];
//        chrono::steady_clock::time_point b4 = chrono::steady_clock::now();
        for (int c = 0; c < l; c++) {
            for (int e = 0; e < counter + 30; e++) {
                if (a[c][e] >= threshold) {
                    int q_number = index_to_query[e];
                    auto pred = [q_number](const query &item) {
                        return item.number == q_number;
                    };
                    auto q = find_if(begin(qq), end(qq), pred);
                    if (q->find_loc == 0) {
                        int query_index = e - q->start_pont;
                        int s = ref_lens[ref_id] + 1;
                        int ref_start_point = c;
                        vector<double> q_dist = query_cmaps[q_number];
                        int l_to_seed = q_dist[query_index - 1] - q_dist[0];
////////////////////////////////////////////////////////
////////////////////////////////// bug Ref start point shoould be changed
//                        int push = 0;
                        while (ref_start_point > 0 &&
                               ref_dist[c - 1] - ref_dist[ref_start_point - 1] <= l_to_seed + band_width) {
                            if (ref_dist[c - 1] - ref_dist[ref_start_point - 1] >= l_to_seed - band_width) {
                                q->seeds_straight[ref_start_point].push_back(make_pair(query_index, c));
//                                push = 1;
                            }
                            ref_start_point--;
                        }
//                        if(push==0){
//                            q->seeds_straight[ref_start_point].push_back(make_pair(query_index, c));
//                        }
                        if (l_to_seed < band_width) {
                            ref_start_point = c + 1;
                            while (ref_start_point < s &&
                                   abs(ref_dist[c - 1] - ref_dist[ref_start_point - 1]) <= l_to_seed + band_width) {
                                q->seeds_straight[ref_start_point].push_back(make_pair(query_index, c));        
                                ref_start_point++;
                            }

                        }

                        l_to_seed = q_dist[query_index + w - 1 - 1] - q_dist[0];
                        ref_start_point = c;
//                        push = 0;
                        if (ref_start_point > 0 && ref_start_point < s) {
                            while (ref_start_point < s &&
                                   ref_dist[ref_start_point - 1] - ref_dist[c - 1] <= l_to_seed + band_width) {
                                if (ref_dist[ref_start_point - 1] - ref_dist[c - 1] >= l_to_seed - band_width) {
                                    q->seeds_reverse[ref_start_point].push_back(make_pair(query_index + w - 1, c));
//                                    push = 1;
                                }
                                ref_start_point++;
                            }
//                            if ( push ==0){
//                                q->seeds_reverse[ref_start_point].push_back(make_pair(query_index + w - 1, c));
//                            }
                            if (l_to_seed < band_width) {
                                ref_start_point = c - 1;
                                while (ref_start_point > 0 &&
                                       abs(ref_dist[ref_start_point - 1] - ref_dist[c - 1]) <= l_to_seed + band_width) {
                                    q->seeds_reverse[ref_start_point].push_back(make_pair(query_index + w - 1, c));
                                    ref_start_point--;
                                }
                            }
                        }
                    }

                }
            }
        }
        destroyTwoDimenArrayOnHeapUsingFree(a, 55000, 8500);
//        chrono::steady_clock::time_point b5 = chrono::steady_clock::now();
//        double b6 = chrono::duration_cast<chrono::milliseconds>(b5 - b4).count() / 1000.;
        // printf("making band: %.3f\n", b6);
        //Inja thread konim
//        if (SV_detection == 1) {
//
//            calculate_seeds_score_in_band_SV(qq, ref_id, ref_cmaps[ref_id], query_cmaps);
//            int ali = 1;
//        } else {
        calculate_seeds_score_in_band(qq, ref_id, ref_cmaps[ref_id], query_cmaps);
//        }
//        chrono::steady_clock::time_point b7 = chrono::steady_clock::now();
//        double b8 = chrono::duration_cast<chrono::milliseconds>(b7 - b5).count() / 1000.;
        // printf("Doing band alignment: %.3f\n", b8);
    }
    // print result;
//    if (SV_detection == 1) {
//        for (auto & q : qq) {
//            vector<double> q_dist = query_cmaps[q.number];
//            map<int, vector<breakpoint>>::reverse_iterator it;
//            float first_score = 1;
//            for (it = q.bp_from_start.rbegin(); it != q.bp_from_start.rend(); it++) {
//                if (it != q.bp_from_start.rbegin()) {
//                    float min_in_index = 1;
//                    vector<breakpoint>::iterator v_t = it->second.begin();
//                    while (v_t != it->second.end()) {
//                        if (v_t->score >= first_score) {
//                            v_t = it->second.erase(v_t);
//                        } else {
//                            if (v_t->score < min_in_index) {
//                                min_in_index = v_t->score;
//                            }
//                            ++v_t;
//                        }
//                    }
//                    if (min_in_index < first_score) {
//                        first_score = min_in_index;
//                    }
//                } else {
//                    for (auto & e : it->second) {
//                        if (e.score < first_score) {
//                            first_score = e.score;
//                        }
//                    }
//                }
//            }
//            first_score = 1.0;
//            map<int, vector<breakpoint>>::iterator it2;
//            for (it2 = q.bp_to_end.begin(); it2 != q.bp_to_end.end(); it2++) {
//                if (it2 != q.bp_to_end.begin()) {
//                    float min_in_index = 1;
//                    vector<breakpoint>::iterator v_t = it2->second.begin();
//                    while (v_t != it2->second.end()) {
//                        if (v_t->score >= first_score) {
//                            v_t = it2->second.erase(v_t);
//                        } else {
//                            if (v_t->score < min_in_index) {
//                                min_in_index = v_t->score;
//                            }
//                            ++v_t;
//                        }
//                    }
//                    if (min_in_index < first_score) {
//                        first_score = min_in_index;
//                    }
//                } else {
//                    for (auto & e : it2->second) {
//                        if (e.score < first_score) {
//                            first_score = e.score;
//                        }
//                    }
//                }
//            }
//            for (auto &pair_start:q.bp_from_start) {
//                for (auto &pair_end:q.bp_to_end) {
//                    if (pair_end.first - pair_start.first == 1 || (pair_end.first > pair_start.first &&
//                                                                   q_dist[pair_end.first - 1] -
//                                                                   q_dist[pair_start.first - 1] <
//                                                                   30000) || (pair_end.first < pair_start.first &&
//                                                                              q_dist[pair_start.first - 1] -
//                                                                              q_dist[pair_end.first - 1] <
//                                                                              20000)) {// ya xeili nazdik be ham hastan
//                        for (auto &start:pair_start.second) {
//                            for (auto &end : pair_end.second) {
//                                if ((start.ref_contig != end.ref_contig ||
//                                     abs(long(start.pos) - long(end.pos)) > 30000) && end.score + start.score < 1.2) {
//                                    if (!((start.ref_contig == 6 && end.ref_contig == 15) ||
//                                          (start.ref_contig == 15 && end.ref_contig == 6)) &&
//                                        !((start.ref_contig == 24 && end.ref_contig == 25) ||
//                                          (start.ref_contig == 25 && end.ref_contig == 24))) {
//                                            svout[number] << q.number << "\t" << start.ref_contig << "\t" << start.dir
//                                                      << "\t"
//                                                      << long(start.pos) << "\t" << start.score << "\t"
//                                                      << pair_start.first << "\t" << end.ref_contig << "\t" << end.dir
//                                                      << "\t"
//                                                      << long(end.pos) << "\t" << end.score << "\t" << pair_end.first
//                                                      << "\n";
//                                        // if (pair_end.first - pair_start.first > 1) {
//                                        //     int diff = pair_end.first - pair_start.first - 1;
//                                        //     if (start.dir == '+') {
//                                        //         svout[number] << q.number << "\t" << start.ref_contig << "\t"
//                                        //                       << start.dir << "\t"
//                                        //                       << long(ref_genome_raw[start.ref_contig][max(0,start.label_pos +
//                                        //                                                                diff - 1)])
//                                        //                       << "\t"
//                                        //                       << start.score << "\t" << pair_start.first + diff << "\t"
//                                        //                       << end.ref_contig << "\t" << end.dir << "\t"
//                                        //                       << long(end.pos) << "\t" << end.score << "\t"
//                                        //                       << pair_end.first << "\n";
//                                        //     } else {
//                                        //         svout[number] << q.number << "\t" << start.ref_contig << "\t"
//                                        //                       << start.dir << "\t"
//                                        //                       << long(ref_genome_raw[start.ref_contig][max(0,start.label_pos -
//                                        //                                                                diff - 1)])
//                                        //                       << "\t"
//                                        //                       << start.score << "\t" << pair_start.first + diff << "\t"
//                                        //                       << end.ref_contig << "\t" << end.dir << "\t"
//                                        //                       << long(end.pos) << "\t" << end.score << "\t"
//                                        //                       << pair_end.first << "\n";
//
//                                        //     }
//                                        //     if (end.dir == '+') {
//                                        //         svout[number] << q.number << "\t" << start.ref_contig << "\t"
//                                        //                       << start.dir
//                                        //                       << "\t" << long(start.pos) << "\t" << start.score << "\t"
//                                        //                       << pair_start.first << "\t" << end.ref_contig << "\t"
//                                        //                       << end.dir << "\t"
//                                        //                       << long(ref_genome_raw[end.ref_contig][max(0,end.label_pos -
//                                        //                                                              diff - 1)]) << "\t"
//                                        //                       << end.score
//                                        //                       << "\t" << pair_end.first - diff << "\n";
//
//                                        //     } else {
//                                        //         svout[number] << q.number << "\t" << start.ref_contig << "\t"
//                                        //                       << start.dir
//                                        //                       << "\t" << long(start.pos) << "\t" << start.score << "\t"
//                                        //                       << pair_start.first << "\t" << end.ref_contig << "\t"
//                                        //                       << end.dir << "\t"
//                                        //                       << long(ref_genome_raw[end.ref_contig][max(0,end.label_pos +
//                                        //                                                              diff - 1)]) << "\t"
//                                        //                       << end.score
//                                        //                       << "\t" << pair_end.first - diff << "\n";
//
//                                        //     }
//                                        // }
//                                    }
//                                }
//                            }
//                        }
//
//                    } else if (pair_end.first > pair_start.first) {
//                        break;
//                    }
//                }
//            }
//        }

//    } else {
    for (auto & q : qq) {
        vector<double> q_dist = query_cmaps[q.number];
        double len_q = q_dist.back() * 1.414;
        double mean = 0.0;
        vector<seedData> seeds;
        vector<answer> v;

//        if (q.ans.empty()) {
//            cout << q.number << " empty ans\n";
//        }

        while (!q.ans.empty()) {
            v.push_back(q.ans.top());
            mean = mean + (1 - (q.ans.top().score / len_q));
            q.ans.pop();
        }
        int p = v.size() - 1;
        if (p == 0) {
            if (v[p].dir=='+'){
                seedData newSeed = seedData(v[p].ref_contig, q.number,v[p].pos-1,p+1 );
                seeds.push_back(newSeed);
            }
            else{
                seedData newSeed = seedData(v[p].ref_contig * - 1, q.number, ref_lens[v[p].ref_contig]-v[p].pos + 1,p + 1);
                seeds.push_back(newSeed);
            }

            // sout[number] << q.number << '\t' << v[p].ref_contig << '\t' << v[p].pos << '\t' << v[p].dir << '\t'
            //              << v[p].score << "\t";
            // for (auto &x:v[p].pairs) {
            //     sout[number] << "(" << x.first << "," << x.second << ")|";
            // }
            // sout[number] << endl;
        } else if (p > 0) {
            mean = mean / v.size();
            if ((1 - (v[p].score / len_q)) > 2.2 * mean) {
                for (int z = p; z > max(p - 10, 0); z--) {
                    if (v[z].dir=='+'){
                        seedData newSeed = seedData(v[z].ref_contig, q.number,v[z].pos-1,p-z+1 );
                        seeds.push_back(newSeed);
                    }
                    else{
                        seedData newSeed = seedData(v[z].ref_contig * - 1, q.number, ref_lens[v[z].ref_contig]-v[z].pos + 1,p -z + 1);
                        seeds.push_back(newSeed);
                    }
                    // sout[number] << q.number << '\t' << v[z].ref_contig << '\t' << v[z].pos << '\t' << v[z].dir
                    //              << '\t'
                    //              << v[z].score << "\t";
                    // for (auto &x:v[z].pairs) {
                    //     sout[number] << "(" << x.first << "," << x.second << ")|";
                    // }
                    // sout[number] << endl;
                }
            } else if ((1 - (v[p].score / len_q)) > 1.7 * mean) {
                for (int z = p; z > max(p - 50, 0); z--) {
                    if (v[z].dir=='+'){
                        seedData newSeed = seedData(v[z].ref_contig, q.number,v[z].pos-1,p-z+1 );
                        seeds.push_back(newSeed);
                    }
                    else{
                        seedData newSeed = seedData(v[z].ref_contig * - 1, q.number, ref_lens[v[z].ref_contig]-v[z].pos + 1,p -z + 1);
                        seeds.push_back(newSeed);
                    }
                    // sout[number] << q.number << '\t' << v[z].ref_contig << '\t' << v[z].pos << '\t' << v[z].dir
                    //              << '\t'
                    //              << v[z].score << "\t";
                    // for (auto &x:v[z].pairs) {
                    //     sout[number] << "(" << x.first << "," << x.second << ")|";
                    // }
                    // sout[number] << endl;
                }
            } else if (1 - (v[p].score / len_q) > 1.5 * mean) {
                for (int z = p; z > max(p - 100, 0); z--) {
                    if (v[z].dir=='+'){
                        seedData newSeed = seedData(v[z].ref_contig, q.number,v[z].pos-1,p-z+1 );
                        seeds.push_back(newSeed);
                    }
                    else{
                        seedData newSeed = seedData(v[z].ref_contig * - 1, q.number, ref_lens[v[z].ref_contig]-v[z].pos + 1,p -z + 1);
                        seeds.push_back(newSeed);
                    }
                    // sout[number] << q.number << '\t' << v[z].ref_contig << '\t' << v[z].pos << '\t' << v[z].dir
                    //              << '\t'
                    //              << v[z].score << "\t";
                    // for (auto &x:v[z].pairs) {
                    //     sout[number] << "(" << x.first << "," << x.second << ")|";
                    // }
                    // sout[number] << endl;
                }
            } else {
                int s_counter = 0;
                while (p >= 0) {
                    s_counter++;
                    if (v[p].dir=='+'){
                        seedData newSeed = seedData(v[p].ref_contig, q.number, v[p].pos-1, s_counter);
                        seeds.push_back(newSeed);
                    }
                    else{
                        seedData newSeed = seedData(v[p].ref_contig * - 1, q.number, ref_lens[v[p].ref_contig]-v[p].pos + 1,s_counter);
                        seeds.push_back(newSeed);
                    }

                    // sout[number] << q.number << '\t' << v[p].ref_contig << '\t' << v[p].pos << '\t' << v[p].dir
                    //              << '\t'
                    //              << v[p].score << "\t";
                    // for (auto &x:v[p].pairs) {
                    //     sout[number] << "(" << x.first << "," << x.second << ")|";
                    // }
                    // sout[number] << endl;
                    p--;
                }
            }
        }
        //Cal aligning function
        molSeedMap[q.number] = seeds;
//        cout << q.number << " m_id_OMF " << seeds.size() << " seeds\n";
    }
    return molSeedMap;
}