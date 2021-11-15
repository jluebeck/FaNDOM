import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="Reference directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-a", "--alignment", help="Alignment Xmap file", required=True)
parser.add_argument("-m", "--mol", help="Query Cmap or Bnx file", required=True)
parser.add_argument("-c", "--chrom", help="Chromosome Version", required=True)
parser.add_argument("-g", "--gene", help="Gene Coodinates Dir", required=False)
args = parser.parse_args()

contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}
rev_contig = {}  # exactly reverse of contig_id_map
for chromosome in contig_id_map:
    rev_contig[contig_id_map[chromosome]] = chromosome


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
    rev_contig = {i:i for i in chr_len.keys()}

    
# HG19
else:
    chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0,
               15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0,
               20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0,
               6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}


scale = 50000



def parse_fda(fdaF):
    fda = {}
    with open(fdaF, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                if line.startswith('0'):
                    mold_id = int(line.split('\t')[2])
                elif line.startswith('2'):
                    elements = line.split('\t')[1:]
                    e = {}
                    for i in elements:
                        i = i.split(':')[1]
                        lab = int(i.split(',')[2])
                        score = int(i.strip().split(',')[3][:-1])
                        e[lab] = score
                    fda[mold_id] = e
    return fda


def parse_bnx(bnxF, keep_length=False):
    moleculeD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                    moleculeD[currKey] = {}
                elif line.startswith('1'):
                    z = 1
                    if keep_length:
                        for x in fields[1:]:
                            moleculeD[currKey][z] = float(x)
                            z += 1
                    else:
                        for x in fields[1:-1]:
                            moleculeD[currKey][z] = float(x)
                            z += 1
    return moleculeD


def parse_xmap(xmapf):
    detailFields = ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "Confidence", "QryLen", "RefLen",
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment"]
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}
    return xmapPair


def find_intersection(intervals, new_interval):
    start, end = new_interval
    for (a, b) in intervals:
        if (a < start < b) or (a < end < b) or (a > start and b < end):
            return 1
    return 0


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


def cluster_alignments():
    svs = {}  # dict of chromosome number to dict of position to alignments
    for chromosome in chr_len.keys():
        svs[chromosome] = {}
        svs[chromosome] = defaultdict(lambda: [], svs[chromosome])
    for id in xmap:
        aln_i = xmap[id]
        chrom = contig_id_map[int(aln_i['RefContigID'])]
        ref_start = float(aln_i['RefStartPos'])
        ref_end = float(aln_i['RefEndPos'])
        index_s = int(ref_start / scale)
        index_e = int(ref_end / scale)
        for j in range(index_s, index_e + 1):
            svs[chrom][j].append(aln_i)  # specify this alignment overlap which this cluster region
    return svs


def indel_detection_alinments(ans, alignment_i, region):
    alignment = alignment_i['Alignment']
    q_id = alignment_i['QryContigID']
    ref_id = alignment_i['RefContigID']
    alignment = alignment.split(')')[:-1]
    padding_limit = 10000
    find_del = []
    find_ins = []
    for i in range(1, 4):  # compare at most lengh window 4 started with 2 in alignment pair j
        for j in range(0, len(alignment) - i):
            ref_1 = int(alignment[j].split(',')[0][1:])  # ref label number
            if region * scale < ref[ref_id][ref_1] <= (region + 1) * scale:
                q_1 = int(alignment[j].split(',')[1])  # query label number
                ref_2 = int(alignment[j + i].split(',')[0][1:])  # ref2 label number
                q_2 = int(alignment[j + i].split(',')[1])  # query2 label number
                ref_dist = abs(ref[ref_id][ref_1] - ref[ref_id][ref_2])
                q_dist = abs(mol[q_id][q_1] - mol[q_id][q_2])
                dif = ref_dist - q_dist
                # check if this happen in not first 25 percent start and end of or 10Kbp
                if mol[q_id][min(q_1, q_2)] > min(0.25 * mol[q_id][len(mol[q_id])], padding_limit) and \
                        mol[q_id][max(q_1, q_2)] < max(0.75 * mol[q_id][len(mol[q_id])],
                                                       mol[q_id][len(mol[q_id])] - padding_limit):
                    if dif < -2000:  # insertion with length greater than 2Kbp
                        if not find_intersection(find_ins, (q_1, q_2)):
                            ans['insertion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
                            find_ins.append((q_1, q_2))
                    if dif > 2000:
                        if not find_intersection(find_del, (q_1, q_2)):
                            ans['deletion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
                            find_del.append((q_1, q_2))
    return ans


def write_indel(selected, indel_type, ref_start_end, chromosome, sv_size):
    len_selected_before_set = len(selected)
    selected = list(set(selected))
    return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(indel_type, chromosome,
                                                        ref[str(rev_contig[chromosome])][
                                                            ref_start_end[0]],
                                                        ref[str(rev_contig[chromosome])][
                                                            ref_start_end[1]], float(
            float(sv_size) / len_selected_before_set), ','.join(selected),
                                                        len(set(selected)))


def indel_detect_concordant_alignment(lines_to_write):
    for chromosome in svs:
        region_to_alignments = svs[chromosome]  # b = ditc of cluster region to list of alignments
        for region in region_to_alignments:  # kk = cluster region number
            ans = {'insertion': [], 'deletion': []}
            for alignment_i in region_to_alignments[region]:  # alignment overlap kk
                ans = indel_detection_alinments(ans, alignment_i, region)
            for indel_type in ans:
                if len(ans[indel_type]) > 0:  # chech is there any indel here in this region?
                    ans[indel_type] = sorted(ans[indel_type], key=lambda x: (x[1], x[2]))
                    selected = []  # Selected SVid
                    ref_start_end = (
                        max(1, ans[indel_type][0][1]),
                        min(ans[indel_type][0][2], len(ref[str(rev_contig[chromosome])])))
                    sv_size = 0  # Sum SV size
                    for aln in ans[indel_type]:
                        # because it is sorted
                        if aln[1] >= ref_start_end[0] - 0 and aln[2] <= ref_start_end[1] + 0:
                            ref_start_end = (min(ref_start_end[0], aln[1]), max(ref_start_end[1], aln[2]))
                            selected.append(aln[0])
                            sv_size += aln[5]
                        else:
                            if len(selected) > 0:
                                lines_to_write.append(
                                    write_indel(selected, indel_type, ref_start_end, chromosome, sv_size))
                            ref_start_end = (max(1, aln[1] - 0), min(aln[2] + 0, len(ref[str(rev_contig[chromosome])])))
                            selected = [aln[0]]
                            sv_size = aln[5]
                    if len(selected) > 0:
                        lines_to_write.append(write_indel(selected, indel_type, ref_start_end, chromosome, sv_size))
    return lines_to_write


def parse_xmap_alignment():
    id_to_partial_alignment = {}  # dict molecule ID to partial alignment of this molecule
    # part2 looking for indel in 2 separate partial alignments
    id_to_partial_alignment = defaultdict(lambda: [], id_to_partial_alignment)
    with open(args.alignment, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line2 = line.strip().split('\t')
                id = int(line2[1])
                id_to_partial_alignment[id].append(line2)
    return id_to_partial_alignment


def detect_indel_discordant(id_to_partial_alignment):
    lines_to_write = []
    for mol_id in id_to_partial_alignment.keys():  # k = molecule ID
        for i in range(len(id_to_partial_alignment[mol_id])):
            q = id_to_partial_alignment[mol_id][i]  # partial alignment 1
            align = q[-1]
            chrom = q[2]
            q_start_label = int(align.split(',')[1].split(')')[0])  # molecule start label number
            ref_start_label = int(align.split(',')[0][1:])  # ref start label number
            ref_start_coordinate = q[5]  # reference start coordiante
            q_end_lanel = int(align.split(',')[-1].split(')')[0])  # molecule end label number
            ref_end_label = int(align.split(',')[-2].split('(')[1])  # ref end label number
            ref_end_coordinate = q[6]  # reference end coordiante
            direction = q[7]
            if direction == '+':
                next_mol_label_start = q_start_label - 1  # next molecule label expected in start
                next_mol_label_end = q_end_lanel + 1  # next molecule label expected in end
                next_mol_coordinate_start = mol[str(mol_id)][
                    max(1, next_mol_label_start)]  # next molecule coordinate expected in start
                next_mol_coordinate_end = mol[str(mol_id)][
                    min(next_mol_label_end, len(mol[str(mol_id)]))]  # next molecule coordinate expected in end
            else:
                next_mol_label_start = q_start_label + 1
                next_mol_label_end = q_end_lanel - 1
                next_mol_coordinate_start = mol[str(mol_id)][min(next_mol_label_start, len(mol[str(mol_id)]))]
                next_mol_coordinate_end = mol[str(mol_id)][max(1, next_mol_label_end)]
            for j in range(len(id_to_partial_alignment[mol_id])):
                if j != i:
                    x1 = 0
                    x2 = 0
                    p = id_to_partial_alignment[mol_id][j]  # partial alignment 2
                    if p[7] == direction and p[2] == chrom:
                        # check have overlap with expected one
                        if min(float(p[3]), float(p[4])) <= next_mol_coordinate_start <= max(float(p[3]),
                                                                                             float(p[4])) or min(
                            float(p[4]), float(p[3])) <= next_mol_coordinate_end <= max(float(p[3]), float(p[4])):
                            align_p = p[-1].split('(')[1:]
                            selected_pair = align_p[:3] + align_p[-3:]
                            for pair in selected_pair:
                                if str(next_mol_label_start) == pair.split(',')[1][:-1]:
                                    p_ref_start = ref[str(chrom)][int(pair.split(',')[0])]
                                    # check if distance is greater than 2KBp and not to faar
                                    if abs(float(ref_start_coordinate) - float(p_ref_start)) - abs(
                                            next_mol_coordinate_start - mol[str(mol_id)][q_start_label]) > 2000 and abs(
                                        float(ref_start_coordinate) - float(p_ref_start)) < 300000:
                                        lines_to_write.append(
                                            "deletion\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(contig_id_map[int(chrom)],
                                                                                              p_ref_start,
                                                                                              ref_start_coordinate, abs(
                                                    float(ref_start_coordinate) - float(p_ref_start)) - abs(
                                                    next_mol_coordinate_start - mol[str(mol_id)][q_start_label]),
                                                                                              mol_id,
                                                                                              1))
                                    x1 = 1
                                if str(next_mol_label_end) == pair.split(',')[1][:-1]:
                                    p_ref_end = ref[str(chrom)][int(pair.split(',')[0])]
                                    if abs(float(ref_end_coordinate) - float(p_ref_end)) - abs(
                                            next_mol_coordinate_end - mol[str(mol_id)][q_end_lanel]) > 2000 and abs(
                                        float(ref_end_coordinate) - float(p_ref_end)) < 300000:
                                        lines_to_write.append(
                                            "deletion\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(contig_id_map[int(chrom)],
                                                                                              ref_end_coordinate,
                                                                                              p_ref_end, abs(
                                                    float(ref_end_coordinate) - float(p_ref_end)) - abs(
                                                    next_mol_coordinate_end - mol[str(mol_id)][q_end_lanel]), mol_id,
                                                                                              1))
                                    x2 = 1
                                if x1 + x2 == 2:
                                    break
    return lines_to_write


def cluster_reads():
    with open(args.output, 'w')as file:
        file.write("#Header\tType\tChromosome\tRefStartPos\tRefEndPos\tPredictedSize\tSupportsId\tSupportsCount\tGenes\n")
        id_to_partial_alignment = {}
        id_to_partial_alignment = defaultdict(lambda: defaultdict(list), id_to_partial_alignment)
        for line_t in lines_to_write:
            c = int(line_t.split('\t')[1])
            sv_type = line_t.split('\t')[0]
            id_to_partial_alignment[c][sv_type].append(line_t)
        for mol_id in id_to_partial_alignment:
            for sv_type in id_to_partial_alignment[mol_id]:
                svs = id_to_partial_alignment[mol_id][sv_type]
                if len(svs) > 0:
                    indel_list = []
                    for i in range(len(svs) - 1):
                        aln = svs[i]
                        nex_aln = svs[i + 1]
                        aln_2 = aln.strip().split('\t')
                        nex_aln_2 = nex_aln.strip().split('\t')
                        if abs(float(aln_2[3]) - float(nex_aln_2[2])) < 30000:
                            nex_aln_2[2] = aln_2[2]
                            nex_aln_2[3] = str(max(float(nex_aln_2[3]), float(aln_2[3])))
                            nex_aln_2[5] += ',' + aln_2[5]
                            nex_aln_2[5] = ','.join(list(set(nex_aln_2[5].split(','))))
                            nex_aln_2[6] = str(len(nex_aln_2[5].split(',')))
                            nex_aln_2[4] = str(float(nex_aln_2[4]) + float(aln_2[4]))
                            svs[i + 1] = '\t'.join(nex_aln_2) + '\n'
                        elif (abs(float(aln_2[2]) - float(nex_aln_2[2])) < 30000 and abs(
                                float(aln_2[3]) - float(nex_aln_2[3])) < 30000):
                            nex_aln_2[2] = aln_2[2]
                            nex_aln_2[3] = str(max(float(nex_aln_2[3]), float(aln_2[3])))
                            nex_aln_2[5] += ',' + aln_2[5]
                            nex_aln_2[5] = ','.join(list(set(nex_aln_2[5].split(','))))
                            nex_aln_2[6] = str(len(nex_aln_2[5].split(',')))
                            nex_aln_2[4] = str(max(float(nex_aln_2[4]), float(aln_2[4])))
                            svs[i + 1] = '\t'.join(nex_aln_2) + '\n'
                        else:
                            indel_list.append(svs[i])
                    indel_list.append(svs[-1])
                    for i in indel_list:
                        i = i.strip().split('\t')
                        chrom = i[1]
                        pos_start = float(i[2])
                        pos_end = float(i[3])
                        if str(chrom) == '23':
                            chrom = 'X'
                        if str(chrom) == '24':
                            chrom = 'Y'
                        chrom = 'chr' + str(chrom)
                        gene_interupt = []
                        genes_chr = genes[chrom]
                        for g in genes_chr:
                            if max(pos_start, g[0]) < min(pos_end, g[1]):
                                gene_interupt.append(g[-1])
                        gene_interupt = list(set(gene_interupt))
                        i.append(','.join(gene_interupt))
                        file.write('\t'.join(i) + '\n')


if __name__ == '__main__':
    ref = parse_cmap(args.ref)
    xmap = parse_xmap(args.alignment)
    mol = args.mol
    if mol.endswith('.cmap'):
        mol = parse_cmap(mol)
    elif mol.endswith('.bnx'):
        mol = parse_bnx(mol)
    if args.chrom == 'nh':
      genes = {}
      genes = defaultdict(lambda: [], genes)
    else:
      genes = parse_gene(args.gene)
    svs = cluster_alignments()
    lines_to_write = []
    lines_to_write = indel_detect_concordant_alignment(lines_to_write)
    id_to_partial_alignment = parse_xmap_alignment()
    lines_to_write = lines_to_write + detect_indel_discordant(id_to_partial_alignment)
    lines_to_write = set(lines_to_write)
    lines_to_write = sorted(lines_to_write, key=lambda x: (
        x.strip().split('\t')[0], int(x.strip().split('\t')[1]), float(x.strip().split('\t')[2]),
        float(x.strip().split('\t')[3])))
    cluster_reads()
