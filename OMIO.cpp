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
    int seed_num = 0;
    ifstream infile(fname);
    if (infile.is_open()) {
        string line;
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
    cout << "Read " << seed_num << " seeds\n";
    return mol_seed_data;
}

void write_xmap_alignment(vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref, map<int,vector<double>> &mol_map,
                          const string &outname, const string &argstring, const bool is_mm) {

    ofstream outfile;
    outfile.open(outname);

    //write header
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    outfile << "# hostname=" << hostname << "\n# $" << argstring << "\n";
    outfile << "# XMAP File Version:\t0.2\n";
    outfile << "# Label Channels:\t1\n";
    outfile << "#h\tXmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\n";
    outfile << "#f\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\tstring\tfloat\tfloat\tint\tstring\n";

    //write alignments
    size_t a_pair_ind = 1;
    for (auto &aln_struct: aln_list) {
        if (!is_mm && aln_struct.is_secondary && !aln_struct.is_partial) {
            continue;
        }
        int raw_ref_id = aln_struct.ref_id;
        int ref_id = abs(raw_ref_id);
        int mol_id = aln_struct.mol_id;
        vector<double> ref_posns = cmaps_ref[ref_id];
        vector<double> mol_posns = mol_map[mol_id];
        vector<tuple<int,int,long>> alignment = aln_struct.alignment;
        float mol_start_pos = mol_posns[get<1>(alignment[0])];
        float mol_end_pos = mol_posns[get<1>(alignment.back())];
        string outline = to_string(a_pair_ind) + "\t" + to_string(mol_id) + "\t" + to_string(ref_id) + "\t";
//        outfile << a_pair_ind << "\t" << mol_id << "\t" << ref_id << "\t" << mol_start_pos << "\t" << mol_end_pos << "\t";
        a_pair_ind++;

        float ref_start_pos;
        float ref_end_pos;
        int ref_start_label;
        long score = get<2>(alignment.back());
        char direction;
        if (raw_ref_id < 0) {
            direction = '-';
            ref_start_label = (ref_posns.size() - get<0>(alignment[0])) - 2;
            int ref_end_label = (ref_posns.size() - get<0>(alignment.back())) - 2;
            ref_start_pos = ref_posns[ref_start_label];
            ref_end_pos = ref_posns[ref_end_label];

        } else {
            direction = '+';
            ref_start_label = get<0>(alignment[0]);
            ref_start_pos = ref_posns[get<0>(alignment[0])];
            ref_end_pos = ref_posns[get<0>(alignment.back())];
        }
        if (direction == '+') {
            outline += (to_string(mol_start_pos) + "\t" + to_string(mol_end_pos) + "\t" + to_string(ref_start_pos) + "\t" + to_string(ref_end_pos) + "\t");
        } else {
            outline += (to_string(mol_end_pos) + "\t" + to_string(mol_start_pos) + "\t" + to_string(ref_end_pos) + "\t" + to_string(ref_start_pos) + "\t");
        }
        outline+=direction;
        outline+=("\t" + to_string(score) + "\t");
        vector<string> alnStringVect;
        vector<string> hitEnum;
        int prev_ref_lab = ref_start_label;
        int prev_mol_lab = get<1>(alignment[0]);
        int m = -1;
        for (const auto &aln_tup: alignment) {
            int d, i = 0;
            m+=1;
            int ref_lab = get<0>(aln_tup);
            int mol_lab = get<1>(aln_tup);
            //translate the reverse label number
            if (raw_ref_id < 0) {
                int n_labs = cmaps_ref[raw_ref_id].size();
                ref_lab = (n_labs - ref_lab) - 2;
            }
            //Build alignment string
            alnStringVect.push_back("(" + to_string(ref_lab+1) + "," + to_string(mol_lab+1) + ")");

            //Build hitEnum
            d=(abs(ref_lab - prev_ref_lab) - 1);
            i=(mol_lab - prev_mol_lab - 1);

            prev_ref_lab = ref_lab;
            prev_mol_lab = mol_lab;

            if (d > 0 or i > 0) {
                hitEnum.push_back(to_string(m) + "M");
                m = 0;
            }

            if (d > 0) {
                hitEnum.push_back(to_string(d) + "D");
            }

            if (i > 0) {
                hitEnum.push_back(to_string(i) + "I");
            }
        }
        if (m > 0) {
            hitEnum.push_back(to_string(m+1) + "M");
        }

        string alnString, hitEnumString;
        if (direction == '-') {
            reverse(alnStringVect.begin(), alnStringVect.end());
            reverse(hitEnum.begin(), hitEnum.end());
        }

        for (const auto &piece : hitEnum) hitEnumString += piece;
        outline+=(hitEnumString + "\t" + to_string(mol_posns.back()) + "\t" + to_string(ref_posns.back()) + "\t1\t");
        for (const auto &piece : alnStringVect) alnString += piece;
        outline+=(alnString + "\n");
        outfile << outline;

    }
    outfile << flush;
    outfile.close();

}

//TODO: REFACTOR
void write_fda_alignment(vector<Alignment> &aln_list, map<int,vector<double>> &cmaps_ref, map<int,vector<double>> &mol_map,
                     const string &outname, const bool is_mm) {

    ofstream outfile;
    outfile.open(outname);

    outfile << "#0\tref_id\tmol_id\taln_direction\tref_start_pos\tref_end_pos\tmol_start_pos\tmol_end_pos\tmol_length\n";
    outfile << "#1\ttotal_score\tmean_score\tis_multimapped\tis_secondary\taln_seed_num\n";
    outfile << "#2\talignment [aln_index]:(ref_pos, mol_pos, mol_lab, score_delta)\n";
    outfile << "#3\tcigar [aln_index]:(delta_ref, delta_mol, mol_label_diff, delta_difference)\n";

    size_t a_pair_ind = 0;
    for (const auto &aln_struct: aln_list) {
        if (!is_mm && aln_struct.is_secondary && !aln_struct.is_partial) {
            continue;
        }
        int raw_ref_id = aln_struct.ref_id;
        int ref_id = abs(raw_ref_id);
        int mol_id = aln_struct.mol_id;
        vector<double> ref_posns = cmaps_ref[ref_id];
        vector<double> mol_posns = mol_map[mol_id];
        vector<tuple<int,int,long>> alignment = aln_struct.alignment;

        outfile << "0\t" << ref_id << "\t" << mol_id << "\t";
        //compute ref start and end posns
//        float ref_start_pos = cmaps_ref[raw_ref_id][get<0>(alignment[0])];
//        float ref_end_pos = cmaps_ref[raw_ref_id][get<0>(alignment.back())];
        float ref_start_pos;
        float ref_end_pos;
        long score = get<2>(alignment.back());
        double mean_score = double(score)/(alignment.size()-1);
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
        long prev_score = 0;
        float prev_ref_pos = ref_start_pos;
        float prev_mol_pos = mol_start_pos;
        int prev_mol_lab = get<1>(alignment[0]);
        int aln_pos = 0;
        for (const auto &aln_tup: alignment) {
            int ref_lab = get<0>(aln_tup);
            int mol_lab = get<1>(aln_tup);
            long curr_score = get<2>(aln_tup);
            long score_delta = curr_score - prev_score;
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



