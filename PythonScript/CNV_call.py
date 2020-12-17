import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-q", "--query", help="Query Dir", required=True)
parser.add_argument("-r", "--ref", help="Reference directory", required=True)
args = parser.parse_args()
contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24,23:25}
from collections import defaultdict

ref_cov = {}
ref_cov = defaultdict(lambda: {}, ref_cov)
query_pairs = {}
query_pairs = defaultdict(lambda:[], query_pairs)
with open(args.ref, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line2 = line.strip().split('\t')
            id = int(line2[0])
            ref_cov[contig_id_map[id]][int(line2[3])] = 0

query_cov = {}
query_cov = defaultdict(lambda: {}, query_cov)
with open(args.query, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line2 = line.strip().split('\t')
            id = int(line2[0])
            query_cov[id][int(line2[3])] = float(line2[7])

with open(args.input, 'r') as f:
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
                if (pair_ref,pair_query,chrom) not in query_pairs[id]:
                	ref_cov[chrom][pair_ref] = ref_cov[chrom][pair_ref] + query_cov[id][pair_query]
                	query_pairs[id].append((pair_ref,pair_query,chrom))
with open(args.output, 'w') as file:
    with open(args.ref, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line2 = line.strip().split('\t')
                id = contig_id_map[int(line2[0])]
                line2[0] = str(id)
                line2[7] = str(ref_cov[id][int(line2[3])])
                line = '\t'.join(line2) + '\n'
            file.write(line)
