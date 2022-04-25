import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output Xmap file", required=True)
args = parser.parse_args()


def parse_xmap_entry(line):
    line2 = line.strip().split('\t')
    contig_id = int(line2[1])
    confidence = float(line2[8])
    query_start_pos = float(line2[3])
    query_end_pos = float(line2[4])
    align = line2[13]
    alignment_pair = align.split((')('))
    alignment_pair[0] = alignment_pair[0][1:]
    query_start_index = int(align.split(')')[0].split(',')[1])
    query_end_index = int(align.split(')')[-2].split(',')[1])
    align = align.split(',')
    number_pairs = len(align) - 1
    score = confidence / number_pairs
    contig_alignments[contig_id].append((1, line))


def filter_contig_alignments():
    for k in contig_alignments.keys():
        deleted = []
        for i in range(len(contig_alignments[k])):
            alignment_i = contig_alignments[k][i]
            alignment_i_line = alignment_i[1].split('\t')
            for j in range(len(contig_alignments[k])):
                if i != j:
                    alignment_j = contig_alignments[k][j]
                    alignment_j_line = alignment_j[1].split('\t')
                    if min(float(alignment_j_line[3]), float(alignment_j_line[4])) <= min(float(alignment_i_line[3]),float(alignment_i_line[4])) and max(float(alignment_j_line[3]),float(alignment_j_line[4])) >= max(
                        float(alignment_i_line[3]), float(alignment_i_line[4])):
                        if abs(float(alignment_j_line[3]) - float(alignment_j_line[4])) - abs(
                                float(alignment_i_line[3]) - float(alignment_i_line[4])) > 50000:
                            deleted.append(alignment_i)
        contig_alignments[k] = [x for x in contig_alignments[k] if not x in deleted]


contig_alignments = {}
contig_alignments = defaultdict(lambda: [], contig_alignments)
if __name__ == '__main__':
    counter = 1
    with open(args.output + '.xmap', 'w') as g:
        with open(args.input, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parse_xmap_entry(line)
                else:
                    g.write(line)
        # for k in contig_alignments.keys():
        #     contig_alignments[k] = sorted(contig_alignments[k], reverse=True)[:300]
        filter_contig_alignments()
        for k in contig_alignments.keys():
            for alignment_i in contig_alignments[k]:
                l = alignment_i[1].strip().split('\t')
                l[0] = str(counter)
                l = '\t'.join(i for i in l)+ '\n'
                g.write(l)
                counter+=1

