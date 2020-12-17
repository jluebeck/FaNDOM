import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-l", "--limit", help="Minimum limit for reporting SV", required=True)
parser.add_argument("-c", "--chrom", help="Chromosome Version", required=True)
parser.add_argument("-q", "--query", help="Query Dir", required=True)
parser.add_argument("-r", "--ref", help="Reference directory", required=True)
args = parser.parse_args()
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
 
mol = parse_cmap(args.query)
ref = parse_cmap(args.ref)

from collections import defaultdict
d = {}
d = defaultdict(lambda:[],d)
with open(args.input, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line2 = line.strip().split('\t')
            id = int(line2[1])
            d[id].append(line2)
a = {}
for i in range(1, 25):
    for j in range(i, 25):
        a[(i, j)] = []

contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}
if args.chrom == 'hg38':
#HG38
    chr_len = {1: 248956422.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0, 15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 2: 242193529.0, 20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0, 6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 23: 156040895.0, 24: 57227415.0}
#HG19
else:
    chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0, 15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0, 20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0, 6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}
limit = int(args.limit)
duplications =[]
inversions = []
for k in d.keys():
    if len(d[k])<200:
        for i in range(len(d[k])):
            q = d[k][i]
            align = q[-1]
            s = (int(align.split(',')[1].split(')')[0]),int(float(q[5])),int(align.split(',')[0][1:]))
            e = (int(align.split(',')[-1].split(')')[0]),int(float(q[6])),int(align.split(',')[-2].split('(')[1]))
            s, e = min(s, e), max(s, e)
            q_chrom = int(q[2])
            q_ref_start = min(float(q[5]),float(q[6]))
            q_ref_end = max(float(q[5]),float(q[6]))
            q_query_start = min(float(q[3]),float(q[4]))
            q_query_end = max(float(q[3]),float(q[4]))
            q_align = q[-1].split(')')[:-1]
            q_align_rev = { int(x.split(',')[1]):int(x.split(',')[0][1:])for x in q_align}
            q_align = {int(x.split(',')[0][1:]) : x.split(',')[1] for x in q_align}
            for j in range(i + 1, len(d[k])):
                p = d[k][j]
                align = p[-1]
                s2 = (int(align.split(',')[1].split(')')[0]),int(float(p[5])),int(align.split(',')[0][1:]))
                e2 = (int(align.split(',')[-1].split(')')[0]),int(float(p[6])),int(align.split(',')[-2].split('(')[1]))
                s2, e2 = min(s2, e2), max(s2, e2)
                p_chrom = int(p[2])
                p_ref_start = min(float(p[5]),float(p[6]))
                p_ref_end = max(float(p[5]),float(p[6]))
                p_query_start = min(float(p[3]),float(p[4]))
                p_query_end = max(float(p[3]),float(p[4]))
                p_align = p[-1].split(')')[:-1]
                p_align_rev = { int(x.split(',')[1]):int(x.split(',')[0][1:])for x in p_align}
                p_align = {int(x.split(',')[0][1:]) : x.split(',')[1] for x in p_align}
                detected = 0
                # print(mol[str(k)][1])
                # if -3 < s2[0] - e[0] < 5 :
                #=============  DUPLICATION FINDER    =========================================
                if p_chrom == q_chrom:# and p[7] == q[7]:
                    if min(p_ref_end,q_ref_end) > max(p_ref_start,q_ref_start) : # this two alignments have overlap
                        ref_intersection =list(set(p_align.keys()) & set(q_align.keys()))
                        duplicated= []
                        for ref_index in ref_intersection:
                            if p_align[ref_index] != q_align[ref_index] :
                                duplicated.append(ref_index)
                        # print(duplicated,'dup')
                        if len(duplicated) > 3:
                            # print(p_chrom,duplicated,k)
                            detected = 1
                            if p[7] == q[7]:
                                duplications.append([str(contig_id_map[p_chrom]),ref[str(p_chrom)][min(duplicated)],ref[str(p_chrom)][max(duplicated)],k,'duplication'])
                            else:
                                duplications.append([str(contig_id_map[p_chrom]),ref[str(p_chrom)][min(duplicated)],ref[str(p_chrom)][max(duplicated)],k,'duplication_inverted'])
                #==============================================================================
                #==========================        INVERSION FINDER       =====================
                if p_chrom == q_chrom and p[7] != q[7] and detected == 0:
                    if max(q_ref_start,p_ref_start) < min(q_ref_end,p_ref_end) or max(q_ref_start,p_ref_start) - min(q_ref_end,p_ref_end) < 100000:
                        if max(q_query_start,p_query_start) < min(q_query_end,p_query_end) or max(q_query_start,p_query_start) - min(q_query_end,p_query_end) < 100000:
                            if p_ref_end < q_ref_end:
                                next_in_ref = max(p_align.keys()) + 1 
                                if next_in_ref in q_align.keys() or ref[str(p_chrom)][min(q_align.keys())] - ref[str(p_chrom)][max(p_align.keys())] < 100000:
                                    inversion_start = ref[str(p_chrom)][max(next_in_ref,min(q_align.keys()))]
                                    last_ind = int(p_align[max(p_align.keys())])+1
                                    if last_ind in q_align_rev.keys():
                                        inversion_end = ref[str(p_chrom)][q_align_rev[last_ind]]
                                    else:
                                        inversion_end = ref[str(p_chrom)][q_align_rev[min(q_align_rev.keys())]]
                                    inversions.append([str(contig_id_map[p_chrom]), inversion_start, inversion_end,k,'inversion'])
                                    detected = 1
                            if q_ref_end < p_ref_end:
                                next_in_ref = max(q_align.keys()) + 1 
                                if next_in_ref in p_align.keys() or ref[str(p_chrom)][min(p_align.keys())] - ref[str(p_chrom)][max(q_align.keys())] < 100000:
                                    inversion_start = ref[str(p_chrom)][max(next_in_ref,min(q_align.keys()))]
                                    last_ind = int(q_align[max(q_align.keys())])+1
                                    if last_ind in p_align_rev.keys():
                                        inversion_end = ref[str(p_chrom)][p_align_rev[last_ind]]
                                    else:
                                        inversion_end = ref[str(p_chrom)][p_align_rev[min(p_align_rev.keys())]]
                                    inversions.append([str(contig_id_map[p_chrom]), inversion_start, inversion_end,k,'inversion'])
                                    detected = 1

                #==============================================================================
                if abs(mol[str(k)][s2[0]] - mol[str(k)][e[0]]) < 30000 or -3 < s2[0] - e[0] < 6:
                        chr1 = (contig_id_map[int(q[2])], e[1])
                        chr2 = (contig_id_map[int(p[2])], s2[1])
                        if chr1[0] < chr2[0] or (chr1[0] == chr2[0] and min(chr2[1], chr1[1]) == chr1[1]):
                            if chr1[1] == int(float(q[6])):
                                orientation_q = '+'
                            else:
                                orientation_q = '-'
                            if chr2[1] == int(float(p[6])):
                                orientation_p = '-'
                            else:
                                orientation_p = '+'
                            a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k,orientation_q,orientation_p))                       
                        else:
                            if chr1[1] == int(float(q[6])):
                                orientation_q = '-'
                            else:
                                orientation_q = '+'
                            if chr2[1] == int(float(p[6])):
                                orientation_p = '+'
                            else:
                                orientation_p = '-'
                            a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k,orientation_p,orientation_q))
                elif s2[0] - e[0] > 0:
                    if abs(float(p[3])- float(p[4]))> 60000 and  abs(float(q[3])- float(q[4]))> 60000 :
                        find = 0
                        for ss in range(len(d[k])):
                            if d[k][ss]!=p and d[k][ss]!=q:
                                if max(float(q[3]),float(q[4])) <float(d[k][ss][3]) < min(float(p[3]),float(p[4])) or max(float(q[3]),float(q[4])) <float(d[k][ss][4]) < min(float(p[3]),float(p[4])):
                                    find = 1
                                    break
                        if find == 0:
                            chr1 = (contig_id_map[int(q[2])], e[1])
                            chr2 = (contig_id_map[int(p[2])], s2[1])
                            if chr1[0] < chr2[0] or (chr1[0] == chr2[0] and min(chr2[1], chr1[1]) == chr1[1]):
                                if chr1[1] == int(float(q[6])):
                                    orientation_q = '+'
                                else:
                                    orientation_q = '-'
                                if chr2[1] == int(float(p[6])):
                                    orientation_p = '-'
                                else:
                                    orientation_p = '+'
                                a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k,orientation_q,orientation_p))                       
                            else:
                                if chr1[1] == int(float(q[6])):
                                    orientation_q = '-'
                                else:
                                    orientation_q = '+'
                                if chr2[1] == int(float(p[6])):
                                    orientation_p = '+'
                                else:
                                    orientation_p = '-'
                                a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k,orientation_p,orientation_q))
                    
                # elif -3 < s[0] - e2[0] < 5 :
                elif abs(mol[str(k)][s[0]] - mol[str(k)][e2[0]]) < 30000 or -3 < s[0] - e2[0] < 6:
                        chr1 = (contig_id_map[int(q[2])], s[1])
                        chr2 = (contig_id_map[int(p[2])], e2[1])
                        if chr1[0] < chr2[0] or (chr1[0] == chr2[0] and min(chr2[1], chr1[1]) == chr1[1]):
                            if chr1[1] == int(float(q[6])):
                                orientation_q = '+'
                            else:
                                orientation_q = '-'
                            if chr2[1] == int(float(p[6])):
                                orientation_p = '-'
                            else:
                                orientation_p = '+'
                            a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k,orientation_q,orientation_p))
                        #     a[(chr2[0], chr1[0])].append((min(chr2[1], chr1[1]),max(chr2[1], chr1[1]),k,orientation_p==orientation_q))
                        else:
                            if chr1[1] == int(float(q[6])):
                                orientation_q = '-'
                            else:
                                orientation_q = '+'
                            if chr2[1] == int(float(p[6])):
                                orientation_p = '+'
                            else:
                                orientation_p = '-'
                            a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k,orientation_p,orientation_q))
                elif s[0] - e2[0] > 0:
                    if abs(float(p[3])- float(p[4]))> 60000 and  abs(float(q[3])- float(q[4]))> 60000 :
                        find = 0
                        for ss in range(len(d[k])):
                            if d[k][ss]!=p and d[k][ss]!=q:
                                if max(float(p[3]),float(p[4])) <float(d[k][ss][3]) < min(float(q[3]),float(q[4])) or max(float(p[3]),float(p[4])) <float(d[k][ss][4]) < min(float(q[3]),float(q[4])):
                                    find = 1
                                    break
                        if find == 0:
                            chr1 = (contig_id_map[int(q[2])], s[1])
                            chr2 = (contig_id_map[int(p[2])], e2[1])
                            if chr1[0] < chr2[0] or (chr1[0] == chr2[0] and min(chr2[1], chr1[1]) == chr1[1]):
                                if chr1[1] == int(float(q[6])):
                                    orientation_q = '+'
                                else:
                                    orientation_q = '-'
                                if chr2[1] == int(float(p[6])):
                                    orientation_p = '-'
                                else:
                                    orientation_p = '+'
                                a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k,orientation_q,orientation_p))
                            #     a[(chr2[0], chr1[0])].append((min(chr2[1], chr1[1]),max(chr2[1], chr1[1]),k,orientation_p==orientation_q))
                            else:
                                if chr1[1] == int(float(q[6])):
                                    orientation_q = '-'
                                else:
                                    orientation_q = '+'
                                if chr2[1] == int(float(p[6])):
                                    orientation_p = '+'
                                else:
                                    orientation_p = '-'
                                a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k,orientation_p,orientation_q))
import numpy as np
scale = 30000
with open(args.output,'w') as file:
    file.write('#Header\tChrom1\tRefPos1\tDirection1\tChrom2\tRefPos2\tDirectio2\tNumberofSupport\tIds\tSupports\n')
    for k in a :
        if len(a[k]) > 0:
            chr1 = k[0]
            chr2 = k[1]
            size1 = int(chr_len[chr1]/scale) + 1
            size2 = int(chr_len[chr2] / scale) + 1
            count_matrix = np.zeros((size1,size2))
            selected = []
            check = {}
            check = defaultdict(lambda:[],check)
            ori = {}
            ori = defaultdict(lambda:{('+','+'):[],('+','-'):[],('-','+'):[],('-','-'):[]},ori)
            for x in a[k]:
                i = int(x[0]/scale)
                j = int(x[1]/scale)
        	    # check[(i,j)].append(x[2])
                ori[(i,j)][(x[3],x[4])].append(x[2])
                count_matrix[i][j]+=1
                if count_matrix[i][j]==limit:
                    selected.append((i,j))

            for ans in selected:
                for kal in ori[(ans[0],ans[1])].keys():
                    if len(set(ori[(ans[0],ans[1])][kal])) > limit-1 :
                        file.write(str(k[0])+'\t'+str(ans[0]*scale+15000)+'\t'+kal[0]+'\t'+str(k[1])+'\t'+str(ans[1]*scale+15000)+'\t'+kal[1]+'\t'+str( count_matrix[ans[0]][ans[1]])+'\t' + str(','.join(str(i) for i in set(ori[(ans[0],ans[1])][kal])))+'\t'+str(len(set(ori[(ans[0],ans[1])][kal])))+'\n')                
    for line in duplications:
        file.write(line[4]+'\t'+str(line[0])+'\t'+str(line[1])+'\t'+str(line[2])+'\t'+str(line[3])+'\n')
    for line in inversions:
        file.write(line[4]+'\t'+str(line[0])+'\t'+str(line[1])+'\t'+str(line[2])+'\t'+str(line[3])+'\n')

     

