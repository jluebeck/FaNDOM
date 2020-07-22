from collections import defaultdict
d = {}
d = defaultdict(lambda:[],d)
with open('assembled_post_process.xmap', 'r') as f:
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

chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0, 15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0, 20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0, 6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}
for k in d.keys():
    for i in range(len(d[k])):
        q = d[k][i]
        align = q[-1]
        s = (int(align.split(',')[1].split(')')[0]),int(float(q[5])),int(align.split(',')[0][1:]))
        e = (int(align.split(',')[-1].split(')')[0]),int(float(q[6])),int(align.split(',')[-2].split('(')[1]))
        s, e = min(s, e), max(s, e)
        for j in range(i + 1, len(d[k])):
            p = d[k][j]
            align = p[-1]
            s2 = (int(align.split(',')[1].split(')')[0]),int(float(p[5])),int(align.split(',')[0][1:]))
            e2 = (int(align.split(',')[-1].split(')')[0]),int(float(p[6])),int(align.split(',')[-2].split('(')[1]))
            s2, e2 = min(s2, e2), max(s2, e2)
            
            if -3 < s2[0] - e[0] < 5 :
                if abs(s2[2]-e[2]) >5 or not q[2]==p[2]:
            # and not(abs(abs(s2[0]-e[0]) - abs(r21 -r12)) < 11 ):
                    # if k == 561 and ( p[0]=='8860' or q[0]=='8860'):
                        # print(s2[2],s2[0],e[2],e[0],'a')
                    chr1 = (contig_id_map[int(q[2])], e[1])
                    chr2 = (contig_id_map[int(p[2])], s2[1])
                    if chr1[0] < chr2[0]:
                        a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k))
                    elif chr1[0] == chr2[0]:
                        a[(chr2[0], chr1[0])].append((min(chr2[1], chr1[1]),max(chr2[1], chr1[1]),k))
                    else:
                        a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k))
            elif -3 < s[0] - e2[0] < 5 :
                # and not(abs(abs(s[0]-e2[0]) - abs(r22 -r11)) < 11):
                if abs(e2[2]-s[2]) > 5 or not q[2]==p[2]:
                    # if k == 561 and ( p[0]=='8860' or q[0]=='8860'):
                        # print(e2[2],s[2],'b')
                    chr1 = (contig_id_map[int(q[2])], s[1])
                    chr2 = (contig_id_map[int(p[2])], e2[1])
                    if chr1[0] < chr2[0]:
                        a[(chr1[0], chr2[0])].append((chr1[1], chr2[1],k))
                    elif chr1[0] == chr2[0]:
                        a[(chr2[0], chr1[0])].append((min(chr2[1], chr1[1]),max(chr2[1], chr1[1]),k))
                    else:
                        a[(chr2[0], chr1[0])].append((chr2[1], chr1[1],k))
import numpy as np
scale = 30000
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
        for x in a[k]:
            i = int(x[0]/scale)
            j = int(x[1]/scale)
	    check[(i,j)].append(x[2])
            count_matrix[i][j]+=1
            if count_matrix[i][j]==1:
                selected.append((i,j))

        for ans in selected:
            print(str(k[0])+'\t'+str(ans[0]*scale+15000)+'\t'+str(k[1])+'\t'+str(ans[1]*scale+15000)+'\t'+str( count_matrix[ans[0]][ans[1]])+'\t' + str(','.join(str(i) for i in check[(ans[0],ans[1])]))+'\t'+str(len(set(check[(ans[0],ans[1])]))))

