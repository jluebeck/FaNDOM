import argparse
from collections import defaultdict
import numpy as np
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-l", "--limit", help="Minimum limit for reporting SV", required=True)
parser.add_argument("-c", "--chrom", help="Chromosome Version", required=True)
parser.add_argument("-q", "--query", help="Query Dir", required=True)
parser.add_argument("-r", "--ref", help="Reference directory", required=True)
parser.add_argument("-g", "--gene", help="Gene Coodinates Dir", required=False)
args = parser.parse_args()
contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}
limit = int(args.limit)

def parse_cmap(cmapf, keep_length=True):
    cmaps = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]
            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                if fD["CMapId"] not in cmaps:
                    cmaps[fD["CMapId"]] = {}
                    # contigCovs[fD["CMapId"]] = {}

                # this is not a good way to parse label channel means color channel
                if fD["LabelChannel"] == "1":
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
                    # contigCovs[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Coverage"])

                elif fD["LabelChannel"] == "0" and keep_length:
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
    return cmaps

def calculate_chr(ref_dir):
    ans = {}
    ref = parse_cmap(ref_dir)
    for k in ref:
        ans[int(k)] = float(list(ref[k].values())[-1])
    return ans


##### develope code for it
if args.chrom == 'hg38':
    # HG38
    chr_len = {1: 248956422.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0,
               15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 2: 242193529.0,
               20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0,
               6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 23: 156040895.0, 24: 57227415.0}

elif args.chrom == 'nh':
    chr_len = calculate_chr(args.ref)
    contig_id_map = {i:i for i in chr_len.keys()}
# HG19
else:
    chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0,
               15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0,
               20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0,
               6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}


class Alignment:
    align = ''
    s = 0
    e = 0
    chrom = 0
    ref_start = 0
    ref_end = 0
    query_start = 0
    query_end = 0
    align_rev = ''
    align_start_query = 0
    direction = ''
    align_end_ref = 0
    id = 0

def parse_xmap(xmap_dir):  # generate dictionary from each molecile id to all alignments
    d = {}
    d = defaultdict(lambda: [], d)
    with open(xmap_dir, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line2 = line.strip().split('\t')
                id = int(line2[1])
                alignment = Alignment()
                q = line2
                align = q[-1]
                s = (int(align.split(',')[1].split(')')[0]), int(float(q[5])), int(align.split(',')[0][1:]))
                e = (
                    int(align.split(',')[-1].split(')')[0]), int(float(q[6])),
                    int(align.split(',')[-2].split('(')[1]))
                alignment.s, alignment.e = min(s, e), max(s, e)
                alignment.chrom = int(q[2])
                alignment.ref_start = min(float(q[5]), float(q[6]))
                alignment.ref_end = max(float(q[5]), float(q[6]))
                alignment.query_start = min(float(q[3]), float(q[4]))
                alignment.query_end = max(float(q[3]), float(q[4]))
                q_align = q[-1].split(')')[:-1]
                alignment.align_rev = {int(x.split(',')[1]): int(x.split(',')[0][1:]) for x in q_align}
                alignment.align = {int(x.split(',')[0][1:]): x.split(',')[1] for x in q_align}
                alignment.align_start_query = float(q[3])
                alignment.direction = q[7]
                alignment.align_end_ref = int(float(q[6]))
                alignment.id = str(id)
                # if id == 2551:
                d[id].append(alignment)
    return d


def duplication_finder(q, p, k):
    detected = False
    if q.chrom == p.chrom:
        if min(p.ref_end, q.ref_end) > max(p.ref_start, q.ref_start):  # this two alignments have overlap
            ref_intersection = list(set(p.align.keys()) & set(q.align.keys()))
            duplicated = []
            for ref_index in ref_intersection:
                if p.align[ref_index] != q.align[ref_index]:
                    duplicated.append(ref_index)
            if len(duplicated) > 3:
                if min(p.query_start, q.query_start) == min(p.align_start_query, q.align_start_query):
                    direction = '+'
                    opp_direction = '-'
                else:
                    direction = '-'
                    opp_direction = '+'
                detected = True
                if p.direction == q.direction:
                    duplications.append([str(contig_id_map[p.chrom]), ref[str(p.chrom)][min(duplicated)],
                                         ref[str(p.chrom)][max(duplicated)], k, 'duplication', direction,
                                         direction])
                else:
                    duplications.append([str(contig_id_map[p.chrom]), ref[str(p.chrom)][min(duplicated)],
                                         ref[str(p.chrom)][max(duplicated)], k, 'duplication_inverted',
                                         direction, opp_direction])
    return detected


def inversion_finder(q, p, k):
    detected = False
    max_gap_between_alignment = 300000
    if p.chrom == q.chrom and p.direction != q.direction:
        if max(q.ref_start, p.ref_start) < min(q.ref_end, p.ref_end) or max(q.ref_start, p.ref_start) - min(
                q.ref_end, p.ref_end) < max_gap_between_alignment:
            if max(q.query_start, p.query_start) < min(q.query_end, p.query_end) or max(q.query_start,
                                                                                        p.query_start) - min(
                q.query_end, p.query_end) < max_gap_between_alignment:
                if min(p.query_start, q.query_start) == min(p.align_start_query, q.align_start_query):
                    direction = '+'
                    opp_direction = '-'
                else:
                    direction = '-'
                    opp_direction = '+'
                detected = True
                if p.ref_end < q.ref_end:
                    next_in_ref = max(p.align.keys()) + 1
                    if next_in_ref in q.align.keys() or ref[str(p.chrom)][min(q.align.keys())] - \
                            ref[str(p.chrom)][max(p.align.keys())] < max_gap_between_alignment:
                        inversion_start = ref[str(p.chrom)][max(next_in_ref, min(q.align.keys()))]
                        last_ind = int(p.align[max(p.align.keys())]) + 1
                        if last_ind in q.align_rev.keys():
                            inversion_end = ref[str(p.chrom)][q.align_rev[last_ind]]
                        else:
                            inversion_end = ref[str(p.chrom)][q.align_rev[min(q.align_rev.keys())]]
                        inversions.append(
                            [str(contig_id_map[p.chrom]), inversion_start, inversion_end, k, 'inversion',
                             direction, opp_direction])
                if q.ref_end < p.ref_end:
                    next_in_ref = max(q.align.keys()) + 1
                    if next_in_ref in p.align.keys() or ref[str(p.chrom)][min(p.align.keys())] - \
                            ref[str(p.chrom)][max(q.align.keys())] < max_gap_between_alignment:
                        inversion_start = ref[str(p.chrom)][max(next_in_ref, min(q.align.keys()))]
                        last_ind = int(q.align[max(q.align.keys())]) + 1
                        if last_ind in p.align_rev.keys():
                            inversion_end = ref[str(p.chrom)][p.align_rev[last_ind]]
                        else:
                            inversion_end = ref[str(p.chrom)][p.align_rev[min(p.align_rev.keys())]]
                        inversions.append(
                            [str(contig_id_map[p.chrom]), inversion_start, inversion_end, k, 'inversion',
                             direction, opp_direction])
    return detected


def find_between(q, p, mol_id):
    distance_between = abs(q.query_end - p.query_start)
    minimum_between_length = 50000
    for i in range(len(d[mol_id])):
        if d[mol_id][i] != p and d[mol_id][i] != q:
            if q.query_end < float(d[mol_id][i].query_start) < p.query_start or q.query_end < float(
                    d[mol_id][i].query_end) < p.query_start or (
                    float(d[mol_id][i].query_start) < q.query_end and float(d[mol_id][i].query_end) > p.query_start):
                distance = min(p.query_start, d[mol_id][i].query_end) - max(q.query_end, d[mol_id][i].query_start)
                if distance > minimum_between_length:
                    return True
    return False


def rev_dir(a):
    if a == '+':
        return '-'
    return '+'


def determine_sv_orientation(q, p, sv1, sv2, mol_id):
    if sv1[0] < sv2[0] or (sv1[0] == sv2[0] and min(sv2[1], sv1[1]) == sv1[1]):
        # if sv1[1] == q.align_end_ref:
        #     orientation_q = '+'
        # else:
        #     orientation_q = '-'
        # if sv2[1] == p.align_end_ref:
        #     orientation_p = '-'
        # else:
        #     orientation_p = '+'
        orientation_q = q.direction
        orientation_p = p.direction
        a[(sv1[0], sv2[0])].append((sv1[1], sv2[1], mol_id, orientation_q, orientation_p))
    else:
        # if sv1[1] == q.align_end_ref:
        #     orientation_q = '-'
        # else:
        #     orientation_q = '+'
        # if sv2[1] == p.align_end_ref:
        #     orientation_p = '+'
        # else:
        #     orientation_p = '-'
        orientation_q = rev_dir(q.direction)
        orientation_p = rev_dir(p.direction)
        a[(sv2[0], sv1[0])].append((sv2[1], sv1[1], mol_id, orientation_p, orientation_q))


def sv_finder(q, p, mol_id):
    partial_alignment_minimum_length = 60000
    max_gap_consecutive_alignments = 30000
    max_overlap_limit= 200000
    # at first check if two alignments have overlap in query or not:
    if not min(q.e[0], p.e[0]) > max(q.s[0], p.s[0]):  # no overlap
        indel_check = 0
        if q.s[0] > p.e[0]:
            p, q = q, p
        # print(p.query_start, q.query_start,abs(abs(q.e[1] - p.s[1]) - abs(q.ref_start - p.ref_end)))
        if p.chrom == q.chrom and p.direction == q.direction:
            if p.direction == '+' and p.ref_start > q.ref_end:
                if abs(abs(q.query_end - p.query_start) - abs(q.ref_end - p.ref_start)) < 300000:
                    indel_check = True
                    return
            elif p.direction == '-'  and p.ref_end < q.ref_start:
                if abs(abs(q.query_end - p.query_start) - abs(q.ref_start - p.ref_end)) < 300000:
                    indel_check = True
                    return

        if abs(mol[str(mol_id)][p.s[0]] - mol[str(mol_id)][q.e[0]]) < max_gap_consecutive_alignments or -3 < p.s[0] - \
                q.e[
                    0] < 6:

            sv1 = (contig_id_map[int(q.chrom)], q.e[1])
            sv2 = (contig_id_map[int(p.chrom)], p.s[1])

            determine_sv_orientation(q, p, sv1, sv2, mol_id)
        elif p.s[0] - q.e[0] > 0:
            if abs(p.query_start - p.query_end) > partial_alignment_minimum_length and abs(
                    q.query_start - q.query_end) > partial_alignment_minimum_length:
                if not find_between(q, p, mol_id):
                    sv1 = (contig_id_map[int(q.chrom)], q.e[1])
                    sv2 = (contig_id_map[int(p.chrom)], p.s[1])
                    determine_sv_orientation(q, p, sv1, sv2, mol_id)
    elif not ((q.e[0] > p.e[0] and q.s[0] < p.s[0]) or (q.e[0] < p.e[0] and q.s[0] > p.s[0])):
        if q.s[0] > p.s[0]:  # not whole subset of the other
            p, q = q, p  # first one q and second one p
        if abs(q.query_end - p.query_start) < max_overlap_limit: #overlap not be very long
            if abs(q.query_start - p.query_start) > partial_alignment_minimum_length and abs(
                    q.query_end - p.query_end) > partial_alignment_minimum_length:
                if q.direction == '+':
                    sv1 = (contig_id_map[int(q.chrom)], q.e[1] - abs(q.query_end - p.query_start))
                    q_dist = abs(mol[str(mol_id)][q.e[0] + 1] - mol[str(mol_id)][q.e[0]]) - abs(
                        p.s[1] - abs(q.e[1] - abs(q.query_end - p.query_start)))
                else:
                    sv1 = (contig_id_map[int(q.chrom)], q.e[1] + abs(q.query_end - p.query_start))
                    q_dist = abs(mol[str(mol_id)][q.e[0] + 1] - mol[str(mol_id)][q.e[0]]) - abs(
                        p.s[1] - abs(q.e[1] + abs(q.query_end - p.query_start)))
                sv2 = (contig_id_map[int(p.chrom)], p.s[1])
                if p.chrom == q.chrom and p.direction == q.direction and abs(q_dist) < 300000:
                    return
                determine_sv_orientation(q, p, sv1, sv2, mol_id)
            if abs(q.query_start - p.query_start) > partial_alignment_minimum_length and abs(
                    q.query_end - p.query_end) > partial_alignment_minimum_length:
                sv1 = (contig_id_map[int(q.chrom)], q.e[1])
                if p.direction == '+':
                    sv2 = (contig_id_map[int(p.chrom)], p.s[1] + abs(q.query_end - p.query_start))
                    q_dist = abs(mol[str(mol_id)][p.s[0] - 1] - mol[str(mol_id)][p.s[0]]) - abs(
                        q.e[1] - abs(p.s[1] + abs(q.query_end - p.query_start)))
                else:
                    sv2 = (contig_id_map[int(p.chrom)], p.s[1] - abs(q.query_end - p.query_start))
                    q_dist = abs(mol[str(mol_id)][p.s[0] - 1] - mol[str(mol_id)][p.s[0]]) - abs(
                        q.e[1] - abs(p.s[1] - abs(q.query_end - p.query_start)))
                if p.chrom == q.chrom and p.direction == q.direction and abs(q_dist) < 300000:
                    return
                determine_sv_orientation(q, p, sv1, sv2, mol_id)


def parse_gene(gene_dir):
    genes = {}
    genes = defaultdict(lambda: [], genes)
    with open(gene_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene_name = line[-4]
            chromosome = line[2]
            strand = line[3]
            start = int(line[4])
            end = int(line[5])
            if chromosome.startswith('chr'):
                genes[chromosome].append([start, end, strand, gene_name])
    return genes


def cluster_sv(k):
    chr1, chr2 = k[0], k[1]
    size1, size2 = int(chr_len[chr1] / scale) + 1, int(chr_len[chr2] / scale) + 1
    count_matrix = np.zeros((size1, size2))
    selected = []
    ori = {}
    ori = defaultdict(lambda: {('+', '+'): [], ('+', '-'): [], ('-', '+'): [], ('-', '-'): []}, ori)
    for x in a[k]:
        i = int(x[0] / scale)
        j = int(x[1] / scale)
        ori[(i, j)][(x[3], x[4])].append(x[2])
        count_matrix[i][j] += 1
        if count_matrix[i][j] == limit:
            selected.append((i, j))
    return selected, ori, count_matrix


def convert_chromosome(chrom):
    if str(chrom) == '23':
        chrom = 'X'
    elif str(chrom) == '24':
        chrom = 'Y'
    chrom = 'chr' + str(chrom)
    return chrom


def check_gene_interuption(chromosome, pos, kal):
    gene_interupt_neighbour = 20000
    genes_chr = genes[chromosome]
    gene_interupt = []
    fuse = []
    for g in genes_chr:
        if g[0] - gene_interupt_neighbour <= pos <= g[1] + gene_interupt_neighbour:
            gene_interupt.append(g[-1])
            if kal[0] != g[2]:
                fuse.append(g[-1])
    return gene_interupt, fuse


def write_sv(k, selected):
    text = ''
    for index in selected:
        for kal in ori[(index[0], index[1])].keys():
            if len(set(ori[(index[0], index[1])][kal])) > limit - 1:
                chrom1 = convert_chromosome(k[0])
                chrom2 = convert_chromosome(k[1])
                pos1 = index[0] * scale + 15000
                pos2 = index[1] * scale + 15000
                gene_interupt1, fuse1 = check_gene_interuption(chrom1, pos1, kal)
                gene_interupt2, fuse2 = check_gene_interuption(chrom2, pos2, kal)
                gene_interupt = list(set(gene_interupt1 + gene_interupt2))
                fused = 'False'
                if len(fuse1) > 0 and len(fuse2) > 0:
                    fused = 'Gene_Fusion'
                sv = "{chr1}\t{pos1}\t{dir1}\t{chr2}\t{pos2}\t{dir2}\t{type}\t{ids}\t{numsup}\t{interupt}\t{fuse}\n".format(
                    chr1=k[0], pos1=index[0] * scale + 15000, dir1=kal[0], chr2=k[1], pos2=index[1] * scale + 15000,
                    dir2=kal[1], type='Unknown', ids=','.join(str(i) for i in set(ori[(index[0], index[1])][kal])),
                    numsup=len(set(ori[(index[0], index[1])][kal])), interupt=','.join(gene_interupt), fuse=fused)
                text = text + sv
    return text


def write_sv2(list_sv):
    sv_text = ''
    for line in list_sv:
        chrom1 = convert_chromosome(str(line[0]))
        pos1 = line[1]
        pos2 = line[2]
        gene_interupt1, _ = check_gene_interuption(chrom1, pos1, [0])
        gene_interupt2, _ = check_gene_interuption(chrom1, pos2, [0])
        gene_interupt = list(set(gene_interupt1 + gene_interupt2))
        sv = "{chr1}\t{pos1}\t{dir1}\t{chr2}\t{pos2}\t{dir2}\t{type}\t{ids}\t{numsup}\t{interupt}\t{fuse}\n".format(
            chr1=line[0], pos1=line[1], dir1=line[5], chr2=line[0], pos2=line[2],
            dir2=line[6], type=line[4], ids=line[3],
            numsup=1, interupt=','.join(gene_interupt), fuse='False')
        sv_text = sv_text + sv
    return sv_text

def parse_bnx(bnxF,keep_length = False):
    moleculeD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                elif line.startswith('1'):
                    if keep_length:
                        moleculeD[currKey] = [float(x) for x in fields[1:]]
                    else:
                        moleculeD[currKey] = {i+1: float(fields[1:-1][i]) for i  in range(len(fields[1:-1]))}

    return moleculeD


if __name__ == '__main__':
    if args.query.endswith('cmap'):
        mol = parse_cmap(args.query)
    elif args.query.endswith('bnx'):
        mol = parse_bnx(args.query)
    ref = parse_cmap(args.ref)
    d = parse_xmap(args.input)
    a = {}
    a = defaultdict(lambda: [], a)
    duplications = []
    inversions = []

    for k in d.keys():
        if len(d[k]) < 200:
            for i in range(len(d[k])):
                q = d[k][i]
                for j in range(i + 1, len(d[k])):
                    p = d[k][j]
                    detected = 0
                    # res2 = False
                    res = duplication_finder(q, p, k)
                    if not res:
                        res = inversion_finder(q, p, k)
                    if not res:
                        sv_finder(q, p, k)
    scale = 30000
    if args.chrom == 'nh':
      genes = {}
      genes = defaultdict(lambda: [], genes)
    else:
      genes = parse_gene(args.gene)
    gene_interupt_neighbour = 20000

    duplications.sort()
    inversions.sort()
    duplications = list(duplications for duplications, _ in itertools.groupby(duplications))
    inversions = list(inversions for inversions, _ in itertools.groupby(inversions))

    with open(args.output, 'w') as file:
        file.write(
            '#Header\tChrom1\tRefPos1\tDirection1\tChrom2\tRefPos2\tDirectio2\tType\tIds\tNumSupports\tGeneInterupt\tGeneFusion\n')
        for k in a:
            if len(a[k]) > 0:
                selected, ori, count_matrix = cluster_sv(k)
                file.write(write_sv(k, selected))
        file.write(write_sv2(duplications))
        file.write(write_sv2(inversions))
