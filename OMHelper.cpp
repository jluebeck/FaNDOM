#include "OMHelper.h"

using namespace std;
//using namespace std;

//--------------------------------------------------------
//Seeding  stuff
//Work with seeds and opening alignment windows

//function to merge intervals
vector<pair<int,int>> mergeIntervals(vector<pair<int,int>> ivect) {
    vector<pair<int, int>> s;
    if (!ivect.empty()) {
//        sort(ivect.begin(), ivect.end(), comparePairs);
        sort(ivect.begin(), ivect.end());
//        for (auto iter: ivect) {
//            cout << "(" << iter.first << ", " << iter.second << ") ";
//        }
//        cout << "\n";
        s.push_back(ivect[0]);
        for (int i = 1; i < ivect.size(); i++) {
            pair<int, int> top = s.back();
            if (top.second < ivect[i].first) {
                s.push_back(ivect[i]);
            } else if (top.second < ivect[i].second) {
                top.second = ivect[i].second;
                s.pop_back();
                s.push_back(top);
            }
        }
    }
    return s;
}

//get the alignment regions from the seeds
//to match with a region, molecule must have 2 or more seeds hitting the region
map<int,vector<seedData>> mol_seeds_to_aln_regions(vector<seedData> &mol_seeds, map<int,vector<double>> &ref_cmaps,
        vector<double> &mol_posns, double padding) {

//    map<int,vector<pair<int,int>>> mol_seed_cmap_regions;
    map<int,vector<seedData>> mol_aln_bounds;
    int mol_lab = 0; //Seeds assumed to be aligned to mol pos 0.
    for (auto &y: mol_seeds) {
        vector<double> curr_ref = ref_cmaps[y.ref_id];
        double ref_0_lab_pos = curr_ref[y.ref_0_lab];
        double mol_lab_pos = mol_posns[mol_lab];
        double left_aln_pos, right_aln_pos;

        left_aln_pos = max((ref_0_lab_pos - mol_lab_pos) - padding, curr_ref[0]+1);
        right_aln_pos = min((ref_0_lab_pos + (mol_posns.back() - mol_lab_pos)) + padding,
                             curr_ref.end()[-2]);

        //construct left label bound
        int lb_lab = lower_bound(curr_ref.begin(),curr_ref.end(), left_aln_pos)- curr_ref.begin() - 1;
        //construct right label bound
        int ub_lab = upper_bound(curr_ref.begin(),curr_ref.end(), right_aln_pos) - curr_ref.begin();
        y.ref_aln_lb = lb_lab;
        y.ref_aln_rb = ub_lab;
        mol_aln_bounds[y.ref_id].emplace_back(y);

    //TODO: SORT, FURTHER MERGING AND FILTERING

    }

    return mol_aln_bounds;
}


////OLD VERSION - MAY REPURPOSE LATER

////get the alignment regions from the seeds
////to match with a region, molecule must have 2 or more seeds hitting the region
//map<int,vector<pair<int,int>>> mol_seeds_to_aln_regions(const map<int,vector<pair<int,int>>> &mol_seeds,
//        map<int,vector<double>> &ref_cmaps, vector<double> &mol_posns, int padding) {
//
//    map<int,vector<pair<int,int>>> mol_seed_cmap_regions;
//    map<int,vector<pair<int,int>>> ref_interval_lists;
//    for (const auto &y: mol_seeds) {
//        int ref_id = y.first;
//        vector<double> curr_cmap = ref_cmaps[ref_id];
//        for (const auto &z: y.second) {
//            float condensed_pos = z.first;
//            int mol_lab = z.second;
//            float mol_lab_pos = mol_posns[mol_lab];
//
//            float left_aln_pos;
//            float right_aln_pos;
//
////            if (padding ) {
////TODO: padding is an INT
//            left_aln_pos = fmax((condensed_pos - mol_lab_pos - padding), curr_cmap[0]+1);
//            right_aln_pos = fmin((condensed_pos + (mol_posns.back() - mol_lab_pos) + padding),
//                                 curr_cmap.end()[-2]);
//
////            } else { //NOT ADVISED FOR ALIGNMENT
////                left_aln_pos = fmax(condensed_pos,curr_cmap[0]+1);
////                right_aln_pos = fmin((condensed_pos + (mol_posns[mol_lab+4] - mol_lab_pos)), curr_cmap.back());
////
////            }
//
//            //construct left bound
//            int lb_lab = (lower_bound(curr_cmap.begin(),curr_cmap.end(),left_aln_pos)- curr_cmap.begin()) - 1;
//            //construct right bound
//            int ub_lab = upper_bound(curr_cmap.begin(),curr_cmap.end(),right_aln_pos) - curr_cmap.begin();
//            ref_interval_lists[ref_id].emplace_back(lb_lab,ub_lab);
//        }
//    }
//    //iterate over the interval lists and merge
//
//    for (const auto &j: ref_interval_lists) {
//        int cmap_id = j.first;
//        vector<pair<int,int>> unmerged_intervals = j.second;
//        vector<pair<int,int>> merged_intervals = mergeIntervals(unmerged_intervals);
//
//        for (const auto &elem: unmerged_intervals) {
//            auto iter_of_elem = find(merged_intervals.begin(),merged_intervals.end(),elem);
//            if (iter_of_elem != merged_intervals.end()) {
//                merged_intervals.erase(iter_of_elem);
//            }
//        }
//        mol_seed_cmap_regions[cmap_id] = merged_intervals;
//
//    }
//
//
//    return mol_seed_cmap_regions;
//}



//-------------------------
//Filtering & Bookkeeping

//get molecule IDs where a full alignment was not found
unordered_set<int> get_remap_mol_ids(const vector<Alignment> &aln_list, map<int, vector<double>> &mol_maps,
        const float remap_prop_cut, const double remap_len_cut) {

    unordered_set<int> seen_mols;
    unordered_set<int> fail_mols;
    for (const auto &a: aln_list) {
        if (! a.is_secondary) {
            int m_id = a.mol_id;
            seen_mols.insert(m_id);
            double molLen = mol_maps[m_id].rbegin()[1];
            double aln_p1 = mol_maps[m_id][get<1>(a.alignment[0])];
            double aln_p2 = mol_maps[m_id][get<1>(a.alignment.back())];
            if ((aln_p2 - aln_p1) / molLen < remap_prop_cut) {
                if (max(aln_p1, molLen - aln_p2) > remap_len_cut) {
                    fail_mols.insert(m_id);
                }
            }
        }
    }
    for (const auto &m: mol_maps) {
        if (seen_mols.find(m.first) == seen_mols.end()) {
            fail_mols.insert(m.first);
        }
    }

    return fail_mols;
}

void filter_mols(map<int,vector<double>> &mol_map, int min_map_lab, int min_map_len) {
    set<int> fail_mols;
    for (const auto &pair: mol_map) {
        int key = pair.first;
        vector<double> posns = pair.second;
        if (posns.size()-1 < min_map_lab || posns.back() < min_map_len) { //Assume size on end
            fail_mols.insert(key);
            continue;
        } else {
            for (int i = 0; i < posns.size() - 1; i++) { //Assume size on end
                double dist_diff = posns[i + 1] - posns[i];
                if (dist_diff > 250000) { //TODO: PARAMETERIZE
                    fail_mols.insert(key);
                    break;
                }
            }
        }
    }

    for (int key: fail_mols) {
        mol_map.erase(key);
    }
}

//Construct the set of CMAPs "reverse complement" and add to the CMAP map
//Take a map of CMAP vectors
//Returns a new CMAP map having both forward and backward elements
map<int,vector<double>> make_reverse_cmap(map<int,vector<double>> &cmap_map) {
    map<int,vector<double>> full_cmap_map;
    for (auto &i: cmap_map) {
        vector<double> curr_posns = i.second;
        int map_id = i.first;
        int map_length = curr_posns.size();
        if (map_length == 0) {
            continue;
        }
        //to do: check if length of curr_posns > 1
        int rev_map_id = -1*map_id;
        full_cmap_map[map_id] = curr_posns;
        full_cmap_map[rev_map_id] = vector<double>();
        if (map_length > 1) {
            for (int j = curr_posns.size() - 2; j > -1; --j) {
                full_cmap_map[rev_map_id].push_back(curr_posns.back() - curr_posns[j]);
            }
        }
        full_cmap_map[rev_map_id].push_back(curr_posns.back());
    }
    return full_cmap_map;
}

//Returns a map of OM ID to length in terms of labels
map<int, int> calculate_length(const map<int, vector<double>> &contigs) {
    map<int, int> contigs_to_length;
    for (auto &contig : contigs) {
        contigs_to_length.insert({contig.first, contig.second.size() - 1});
    }
    return contigs_to_length;
}

//-------------------------------------------------------------------------------