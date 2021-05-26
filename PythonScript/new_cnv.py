from collections import defaultdict
import argparse
import numpy as np
import json
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--contig_query", help="contigs cmap file", required=True)
parser.add_argument("-f", "--folderalignment", help="Folder molecule to contigs alignments", required=True)
parser.add_argument("-a", "--contig_alignments", help="contigs alignments", required=True)
parser.add_argument("-r", "--ref", help="reference cmap file", required=True)
parser.add_argument("-o", "--output", help="output directory", required=True)
args = parser.parse_args()
contigs_dir = args.contig_query
alignment_dir = args.folderalignment
contigs_alignment = args.contig_alignments
reference_dir = args.ref
output_dir = args.output

contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24, 23: 25}
contigs = {}
contigs = defaultdict(lambda: {}, contigs)
ids = []
with open(contigs_dir, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            id = int(line[0])
            length = int(line[2])
            ids.append((id, length))

ids = list(set(ids))

for id in ids:
    # if id[0] == 434:
    file_dir = alignment_dir + '/EXP_REFINEFINAL1_noOutlier_contig' + str(id[0]) + '.xmap'
    contigs[id[0]] = defaultdict(lambda: [], contigs[id[0]])
    with open(file_dir, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                molecule_id = int(line[1])
                pair_alignment = line[-1]
                pair_alignment = pair_alignment.split('(')[1:]
                for i in pair_alignment:
                    i = int(i.split(',')[0])
                    contigs[id[0]][i].append(molecule_id)

for key in contigs:
    for key2 in contigs[key]:
        contigs[key][key2] = list(set(contigs[key][key2]))

# with open(alignment_dir+'/data2.json', 'w') as fp:
#     json.dump(contigs, fp)
###########################################################################################
ref_cov = {}
ref_cov = defaultdict(lambda: {}, ref_cov)
query_pairs = {}
query_pairs = defaultdict(lambda: [], query_pairs)

with open(reference_dir, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line2 = line.strip().split('\t')
            id = int(line2[0])
            ref_cov[contig_id_map[id]][int(line2[3])] = []

with open(contigs_alignment, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line2 = line.strip().split('\t')
            id = int(line2[1])
            chrom = contig_id_map[int(line2[2])]
            align = line2[-1]
            align = align.split(')')[:-1]
            for pair in align:
                pair = pair.split(',')
                pair_ref = int(pair[0][1:])
                pair_query = int(pair[1])
                if (pair_ref, pair_query, chrom) not in query_pairs[id]:
                    ref_cov[chrom][pair_ref] = list(set(ref_cov[chrom][pair_ref] + contigs[id][pair_query]))
                    query_pairs[id].append((pair_ref, pair_query, chrom))

# print(ref_cov)
all_cov = []
for chrom in ref_cov.keys():
    for pos in ref_cov[chrom].keys():
        all_cov.append(len(ref_cov[chrom][pos]))
mean = np.mean(all_cov)
std = np.std(all_cov)


with open(output_dir, 'w') as file:
    with open(reference_dir, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line2 = line.strip().split('\t')
                id = contig_id_map[int(line2[0])]
                line2[0] = str(id)
                line2[7] = str(len(ref_cov[id][int(line2[3])]))
                copynumner = str(2 * int(line2[7])/mean)
                if line2[7] == '0':
                    copynumner = '0'
                line2.append(copynumner)
                # line.append(str(2 * (int(line[7])) / mean))
                ref_cov[id][int(line2[3])] = [str(i) for i in ref_cov[id][int(line2[3])]]
                line2.append(','.join(ref_cov[id][int(line2[3])]))
                line = '\t'.join(line2) +'\n'
            file.write(line)
