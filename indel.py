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
                        moleculeD[currKey][z] = x
                        z += 1
                    # moleculeD[currKey] = [float(x) for x in fields[1:]]
                else:
                    for x in fields[1:-1]:
                        moleculeD[currKey][z] = x
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


contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}
rev_contig = {}
for k in contig_id_map:
    rev_contig[contig_id_map[k]] = k

ref = parse_cmap('/home/sraeisid/Merged_FaNDOM/FaNDOM/reference_genomes/hg38_Merge_800_DLE.cmap')
xmap = parse_xmap('/nucleus/projects/sraeisid/OMFilter_result/SV/K562/res/Fandom/sv/hg38/all_hg38.xmap')

mol = '/nucleus/projects/sraeisid/OMFilter_result/SV/K562/molecule_mapping/exp_alignmolvrefsv_q.cmap'

if mol.endswith('.cmap'):
    mol = parse_cmap(mol)
elif mol.endswith('.bnx'):
    mol = parse_bnx(mol)
 #HG38
chr_len = {1: 248956422.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0, 15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 2: 242193529.0, 20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0, 6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 23: 156040895.0, 24: 57227415.0}
# HG19
#chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0,
#           15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0,
 #          20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0,
  #         6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}
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
                    q_dist = abs(mol[q_id][q_1] - mol[q_id][q_2])
                    dif = ref_dist - q_dist
                    if dif < -3000:
                        ans['insertion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
                    if dif > 3000:
                        ans['deletion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])


        for i in ans:
            if len(ans[i]) > 0 :
                ans[i] = sorted(ans[i], key = lambda x : x[1])
                selected = []
                i_j = (ans[i][0][1],ans[i][0][2])
                for aln in ans[i]:
                    if aln[1] == i_j[0] and aln[2] == i_j[1]:
                        selected.append(aln[0])
                    else:
                        if len(selected) > 4:
                            print(i + '\t' + str(k) + '\t'+ str(ref[str(rev_contig[k])][i_j[0]])+ '\t' + str(ref[str(rev_contig[k])][i_j[1]]) + '\t'+','.join(selected)+ '\t' + str(len(selected)))
                        i_j = (aln[1],aln[2])
                        selected = [aln[0]]
                if len(selected) > 4 :
                    print(i + '\t' + str(k) + '\t' + str(ref[str(rev_contig[k])][i_j[0]]) + '\t' + str(
                    ref[str(rev_contig[k])][i_j[1]]) + '\t' + ','.join(selected)+ '\t' + str(len(selected)))


