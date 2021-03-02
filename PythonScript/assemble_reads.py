from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for Xmap file", required=True)
args = parser.parse_args()


def compare_alignment(alignment1, alignment2):
    if alignment2 in alignment1:
        return alignment1
    alignment1 = alignment1.split(')')
    alignment2 = alignment2.split(')')
    while '' in alignment1:
        alignment1.remove('')
    while '' in alignment2:
        alignment2.remove('')
    assembled_alignments = alignment1 + alignment2
    assembled_alignments = list(set(assembled_alignments))
    assembled_alignments = sorted(assembled_alignments, key=lambda x: int(x.split(',')[0][1:]))
    assembled_alignments = ')'.join(assembled_alignments) + ')'
    return assembled_alignments


def join_string(alignment1, alignment2, direction):
    if alignment2 in alignment1:
        return alignment1
    if direction == '+':
        while alignment2 != '':
            alignment_pair = alignment2.split(')')[0]
            alignment_pair = alignment_pair + ')'
            if alignment_pair in alignment1:
                break
            else:
                alignment2 = alignment2.split(')')[1:]
                alignment2 = ')'.join(alignment2)
        while alignment1 != '':
            alignment_pair = alignment1.split(')')[-2]
            alignment_pair = alignment_pair + ')'
            if alignment_pair in alignment2:
                break
            else:
                alignment1 = alignment1[:alignment1.rfind('(')]
        return compare_alignment(alignment1, alignment2)
    else:
        while True:
            alignment_pair = alignment2.split(')')[-2]
            alignment_pair = alignment_pair + ')'
            if alignment_pair in alignment1:
                break
            else:
                alignment2 = alignment2[:alignment2.rfind('(')]
        while True:
            alignment_pair = alignment1.split(')')[0]
            alignment_pair = alignment_pair + ')'
            if alignment_pair in alignment2:
                break
            else:
                alignment1 = alignment1.split(')')[1:]
                alignment1 = ')'.join(alignment1)
        return compare_alignment(alignment1, alignment2)


def parse_input_xmap():
    d = {}
    d = defaultdict(lambda: [], d)
    with open(args.input, 'r')as f:
        for line in f:
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                d[line2[1]].append(line2)
    return d


def remove_subset_alignments(contig_id):
    contigs_alignments[contig_id] = sorted(contigs_alignments[contig_id], key=lambda x: (
        min(float(x[3]), float(x[4])), -max(float(x[3]), float(x[4]))))
    deleted = []
    for i in range(len(contigs_alignments[contig_id])):
        i1 = contigs_alignments[contig_id][i]
        for j in range(len(contigs_alignments[contig_id])):
            if i < j:
                j1 = contigs_alignments[contig_id][j]
                if j1[-1] in i1[-1]:
                    if j1 not in deleted:
                        deleted.append(j1)
    for i in deleted:
        contigs_alignments[contig_id].remove(i)


if __name__ == '__main__':
    with open(args.output + '_assembled.xmap', 'w') as g:
        contigs_alignments = parse_input_xmap()
        for contig_id in contigs_alignments:
            print("assembling contigs {id}".format(id=contig_id))
            remove_subset_alignments(contig_id)
            for i in range(len(contigs_alignments[contig_id])):
                find = 0
                main_aln = contigs_alignments[contig_id][i]
                for j in range(i + 1, len(contigs_alignments[contig_id])):
                    nex_aln = contigs_alignments[contig_id][j]
                    if main_aln[2] == nex_aln[2] and main_aln[7] == nex_aln[7]:  # same chromosome and same direction
                        main_aln_query_end_pos = max(float(main_aln[3]), float(main_aln[4]))
                        nex_aln_query_start_pos = min(float(nex_aln[3]), float(nex_aln[4]))
                        if main_aln_query_end_pos > nex_aln_query_start_pos and len(
                                set(main_aln[-1].split(')')) & set(nex_aln[-1].split(')'))) > 1:
                            if (max(float(main_aln[5]), float(main_aln[6])) > min(float(nex_aln[5]),
                                                                                  float(nex_aln[6])) > min(
                                float(main_aln[5]), float(main_aln[6])) and main_aln[7] == '+') or (
                                    max(float(nex_aln[5]), float(nex_aln[6])) > min(float(main_aln[5]),
                                                                                    float(main_aln[6])) > min(
                                float(nex_aln[5]), float(nex_aln[6])) and main_aln[7] == '-'):
                                if main_aln[7] == '+':
                                    main_aln[4] = nex_aln[4]
                                    main_aln[6] = nex_aln[6]
                                    main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '+')
                                else:
                                    main_aln[3] = nex_aln[3]
                                    main_aln[5] = nex_aln[5]
                                    main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '-')
                                # i += 1
                                find = 1
                                contigs_alignments[contig_id][j] = main_aln
                                break
                        elif main_aln[7] == '+':
                            main_aln_query_end = int(main_aln[-1].split(',')[-1][:-1])
                            main_aln_ref_end = int(main_aln[-1].split(',')[-2].split('(')[1])
                            next_aln_query_start = int(nex_aln[-1].split(',')[1].split(')')[0])
                            next_aln_ref_start = int(nex_aln[-1].split(',')[0][1:])
                            if next_aln_query_start - main_aln_query_end == next_aln_ref_start - main_aln_ref_end == 1:
                                main_aln[4] = nex_aln[4]
                                main_aln[6] = nex_aln[6]
                                main_aln[-1] = main_aln[-1] + nex_aln[-1]
                                find = 1
                                contigs_alignments[contig_id][j] = main_aln
                                break
                                # i+=1
                        elif main_aln[7] == '-':
                            main_aln_query_end = int(main_aln[-1].split(',')[1].split(')')[0])
                            main_aln_ref_end = int(main_aln[-1].split(',')[0][1:])
                            next_aln_query_start = int(nex_aln[-1].split(',')[-1][:-1])
                            next_aln_ref_start = int(nex_aln[-1].split(',')[-2].split('(')[1])
                            if next_aln_query_start - main_aln_query_end == main_aln_ref_end - next_aln_ref_start == 1:
                                main_aln[3] = nex_aln[3]
                                main_aln[5] = nex_aln[5]
                                main_aln[-1] = nex_aln[-1] + main_aln[-1]
                                find = 1
                                contigs_alignments[contig_id][j] = main_aln
                                break
                                # i+=1
                if find == 0:
                    g.write('\t'.join(main_aln))
                    g.write('\n')
                    # main_aln = nex_aln
                    # i += 1
            # g.write('\t'.join(main_aln))
            # g.write('\n')
