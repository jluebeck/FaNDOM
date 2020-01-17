//#include <cstdio>
#include <iostream>
//#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <set>
#include <sstream>
#include "OMIO.h"

using namespace std;


//Parse CMAP file where id is in column 0, and position is in column 5.
//Takes the CMAP filename as input.
//Returns a map of ids to a vector of doubles. The length field IS kept.
map<int,vector<double>> parse_cmap(const string &fname) {
    map<int,vector<double>> cmap_pairs;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
        size_t pos = 0;
        int curr_id = 0;
        double curr_pos = 0;
        while (getline (infile,line)) {
            int counter = 0;
            if (line[0] == '#') {
                continue;
            }
            while ((pos = line.find('\t')) != string::npos) {
                string token = line.substr(0,pos);
                if (counter == 0) {
                    if (curr_id != stoi(token)) {
                        curr_id = stoi(token);
                    }
                } else if (counter == 5) {
                    curr_pos = stof(token);
                    break;
                }
                line.erase(0,pos+1);
                counter++;
            }
            cmap_pairs[curr_id].push_back(curr_pos);
        }
        infile.close();
    }
    //cmap_pairs[curr_id].pop_back(); //REMOVES THE LAST LABEL FROM LAST ENTRY
    return cmap_pairs;
}

//Print out a CMAP map element by element
void cmap_map_to_string(map<int,vector<double>> &cmap_map) {
    for (auto &iter: cmap_map) {
        cout << "CMAP ID " << iter.first << "\n";
        vector<double> curr_posns = iter.second;
        for (auto elem: curr_posns) {
            cout << "POSN " << elem << "\n";
        }
    }
    cout << endl;
}

//map<int,float>parse_score_thresholds(string &score_thresh_file) {
//    map<int,float>score_thresholds;
//    ifstream infile(score_thresh_file);
//    int seg_id;
//    float score;
//    while (infile >> seg_id >> score) {
//        if (score > 0) {
//            score_thresholds[seg_id] = score;
//            score_thresholds[-1 * seg_id] = score;
//        }
//    }
//    return score_thresholds;
//}

set<int>parse_contig_list(string &contig_list_file) {
    set<int>contig_set;
    ifstream infile(contig_list_file);
    int contig_id;
    while (infile >> contig_id) {
        contig_set.insert(contig_id);
    }
    return contig_set;
}

map<int,vector<pair<int,int>>> parse_bed(const string &fname, map<string,pair<int,double>> &key_data) {
    map<int,vector<pair<int,int>>> bed_data;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
        while (getline (infile,line)) {
            istringstream tokenStream(line);
            string token;
            string chr;
            int start,end;
            int counter = 0;
            while (getline (tokenStream,token,'\t')) {
                switch (counter) {
                    case 0:
                        chr = token;
                    case 1:
                        start = stoi(token);
                    case 2:
                        end = stoi(token);
                    default:
                        break;
                }
                counter++;
            }
            if (counter < 3) {
                cout << "Malformatted bed file, exiting.";
                exit(1);
            } else {
                auto val_ptr = key_data.find(chr);
                if (val_ptr != key_data.end()) {
                    int c_id;
                    double len;
                    tie(c_id,len) = val_ptr->second;
                    bed_data[c_id].emplace_back(start,end);

                    //add the reverse
                    bed_data[-c_id].emplace_back(len-end,len-start);

                } else {
                    cout << "Unfound cmap id for bed entry " << line << "\n";
                }
            }
        }
    }
    return bed_data;
}

map<string,pair<int,double>> parse_key(const string &fname) {
    map<string,pair<int,double>> key_data;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
        while (getline (infile,line)) {
            if (line[0] == '#' || line.substr(0,8) == "CompntId" ) {
                continue;
            }
            istringstream tokenStream(line);
            int id;
            string token,chr;
            double len;
            int counter = 0;
            while (getline (tokenStream,token,'\t')) {
                switch (counter) {
                    case 0:
                        id = stoi(token);
                    case 1:
                        chr = token;
                    case 2:
                        len = stod(token);
                    default:
                        break;
                }
                counter++;
            }
            if (counter < 3) {
                cout << "Malformatted key file, exiting.";
                exit(1);
            } else {
                key_data[chr] = {id,len};
            }
        }
    }
    return key_data;
}

map<int,vector<double>> parse_bnx(const string &fname) {
    map<int,vector<double>> bnx_map;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
        size_t pos;
        double curr_pos;
        int curr_id = -1;
        while (getline (infile,line)) {
            if (line[0] == '0') {
                int first_tab = line.find('\t');
                line.erase(0,first_tab+1);
                int second_tab = line.find('\t');
                curr_id = stoi(line.substr(0,second_tab));

            } else if (line[0] == '1') {
                int counter = 0;
                string token;
                while ((pos = line.find('\t')) != string::npos) {
                    token = line.substr(0, pos);
                    if (counter != 0) {
                        curr_pos = stod(token);
                        if (curr_id == -1) {
                            cout << "Error: Malformed BNX." << endl;
                            exit(1);
                        }
                        bnx_map[curr_id].push_back(curr_pos);

                    }
                    line.erase(0,pos+1);
                    counter++;
                }
                //includes the last (size) label
                token = line.substr(0, pos);
                curr_pos = stod(token);
                bnx_map[curr_id].push_back(curr_pos);
            }
        }
    }
    return bnx_map;
}

map<int, vector<seedData>> parse_seeds(const string &fname, map<int, vector<double>> &ref_cmaps) {
    cout << "Parsing seeds from " << fname << "\n";
    map<int, vector<seedData>> mol_seed_data;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
        int seed_num = 0;
        while (getline (infile,line)) {
            if (line[0] == '#') {
                cout << "skipping header line on seeds\n";
                continue;
            }
            seed_num+=1;
            istringstream tokenStream(line);
            int mol_id = 0;
            int ref_id = 0;
            int ref_0_lab = 0;
            string dir = "+";
            string token;
            int counter = 0;
            while (getline (tokenStream,token,'\t')) {
                switch (counter) {
                    case 0:
                        mol_id = stoi(token);
                    case 1:
                        ref_id = stoi(token);
                    case 2:
                        ref_0_lab = stoi(token) - 1; // ASSUMES PROVIDED IN STANDARD CMAP 1-index
                    case 3:
                        dir = token;
                    default:
                        break;
                }
                counter++;
            }
            if (counter < 3) {
                cout << "Malformatted seeds file, exiting.";
                exit(1);

            } else if (mol_id != 0) {
                int ref_len = ref_cmaps[ref_id].size();
                if (dir == "-") {
                    ref_id *= -1;
                    ref_0_lab = (ref_len - ref_0_lab) - 1;
                }

                ref_0_lab = max(0,ref_0_lab);
                ref_0_lab = min(ref_0_lab,ref_len);
                seedData newSeed = seedData(ref_id,mol_id,ref_0_lab,seed_num);
                mol_seed_data[mol_id].push_back(newSeed);
            }
        }
    }
    return mol_seed_data;
}




