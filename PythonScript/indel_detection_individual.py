def parse_cmap(cmapf, keep_length=False):
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
import statistics

def parse_fda(fdaF):
    fda = {}
    with open(fdaF, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                if line.startswith('0'):
                    # print('0',line)
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
                        # moleculeD[currKey] = [float(x) for x in fields[1:]]
                    else:
                        for x in fields[1:-1]:
                            moleculeD[currKey][z] = float(x)
                            z += 1
                        # moleculeD[currKey] = [float(x) for x in fields[1:-1]]
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
                alnstring = ")" + fD["Alignment"] + "("
                # xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}

    return xmapPair



import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="Reference directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-a", "--alignment", help="Alignment Xmap file", required=True)
parser.add_argument("-m", "--mol", help="Query Cmap or Bnx file", required=True)
parser.add_argument("-c", "--chrom", help="Chromosome Version", required=True)
args = parser.parse_args()

contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}
rev_contig = {}
for k in contig_id_map:
    rev_contig[contig_id_map[k]] = k

ref = parse_cmap(args.ref)
xmap = parse_xmap(args.alignment)

mol = args.mol

if mol.endswith('.cmap'):
    mol = parse_cmap(mol)
elif mol.endswith('.bnx'):
    mol = parse_bnx(mol)
# HG38
if args.chrom == '38':
    chr_len = {1: 248956422.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0, 15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 2: 242193529.0, 20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0, 6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 23: 156040895.0, 24: 57227415.0}
# HG19
else:
    chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0,
               15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0,
               20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0,
               6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}
from collections import defaultdict

a = {}

for k in chr_len.keys():
    a[k] = {}
    a[k] = defaultdict(lambda: [], a[k])

scale = 50000
for k in xmap:
    i = xmap[k]
    chrom = contig_id_map[int(i['RefContigID'])]
    ref_start = float(i['RefStartPos'])
    ref_end = float(i['RefEndPos'])
    index_s = int(ref_start / scale)
    index_e = int(ref_end / scale)
    for j in range(index_s, index_e + 1):
        a[chrom][j].append(i)
with open(args.output,'w')as file :
    for k in a:
        b = a[k]
        for kk in b:
            ans = {'insertion': [], 'deletion': []}
            for i in b[kk]:
                alignment = i['Alignment']
                q_id = i['QryContigID']
                ref_id = i['RefContigID']
                alignment = alignment.split(')')[:-1]
                for v in range(0, len(alignment) - 1):
                    ref_1 = int(alignment[v].split(',')[0][1:])
                    if kk * scale < ref[ref_id][ref_1] <= (kk + 1) * scale:
                        q_1 = int(alignment[v].split(',')[1])
                        ref_2 = int(alignment[v + 1].split(',')[0][1:])
                        q_2 = int(alignment[v + 1].split(',')[1])
                        ref_dist = abs(ref[ref_id][ref_1] - ref[ref_id][ref_2])
                        if abs(q_2-q_1) < 3 and abs(ref_2-ref_1)<3:
                            q_dist = abs(mol[q_id][q_1] - mol[q_id][q_2])
                            dif = ref_dist - q_dist
                            if dif < -2000:
                                if mol[q_id][min(q_1, q_2)] > min(0.25 * mol[q_id][len(mol[q_id])],50000) and mol[q_id][max(q_1, q_2)] < max(0.75 * mol[q_id][len(mol[q_id])],  mol[q_id][len(mol[q_id])]-50000):
                                    ans['insertion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
                            if dif > 2000:
                                if mol[q_id][min(q_1, q_2)] > min(0.25 * mol[q_id][len(mol[q_id])],50000) and mol[q_id][max(q_1, q_2)] < max(0.75 * mol[q_id][len(mol[q_id])], mol[q_id][len(mol[q_id])]-50000):
                                    ans['deletion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
            for i in ans:
                if len(ans[i]) > 0:
                    ans[i] = sorted(ans[i], key=lambda x: (x[1], x[2]))
                    selected = []
                    i_j = (max(1,ans[i][0][1]),min(ans[i][0][2],len(ref[str(rev_contig[k])])))
                    s = 0
                    s2 = 0
                    for aln in ans[i]:
                        if aln[1] >= i_j[0] - 3 and aln[2] <= i_j[1] + 3:
                            i_j = (min(i_j[0], aln[1]), max(i_j[1], aln[2]))
                            selected.append(aln[0])
                            s += aln[5]
                            s2 += aln[6]
                        else:
                            if len(selected) > 0:
                                selected = list(set(selected))
                                file.write(i + '\t' + str(k) + '\t' + str(ref[str(rev_contig[k])][i_j[0]]) + '\t' + str(
                                    ref[str(rev_contig[k])][i_j[1]]) + '\t' + ','.join(selected) + '\t' + str(
                                    len(set(selected))) + '\t' + str(float(float(s) / len(selected))) + '\t' + str(
                                    float(float(s2) / len(selected)))+'\n')
                            i_j = (max(1, aln[1] - 3), min(aln[2] + 3, len(ref[str(rev_contig[k])])))
                            selected = [aln[0]]
                            s = aln[5]
                            s2 = aln[6]
                    if len(selected) > 0:
                        selected = list(set(selected))
                        file.write(i + '\t' + str(k) + '\t' + str(ref[str(rev_contig[k])][i_j[0]]) + '\t' + str(
                            ref[str(rev_contig[k])][i_j[1]]) + '\t' + ','.join(selected) + '\t' + str(
                            len(set(selected))) + '\t' + str(float(float(s) / len(selected))) + '\t' + str(
                            float(float(s2) / len(selected)))+'\n')
