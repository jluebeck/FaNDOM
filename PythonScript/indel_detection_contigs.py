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


def find_intersection(intervals, new_interval):
	start, end = new_interval

	for (a, b) in intervals:
		if (a < start < b) or (a < end < b) or (a > start and b < end):
			return 1
	return 0


import argparse

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
rev_contig = {} #exactly reverse of contig_id_map
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
	chr_len = {1: 248956422.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0,
			   15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 2: 242193529.0,
			   20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0,
			   6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 23: 156040895.0, 24: 57227415.0}
# HG19
else:
	chr_len = {1: 249250621.0, 10: 135534747.0, 11: 135006516.0, 12: 133851895.0, 13: 115169878.0, 14: 107349540.0,
			   15: 102531392.0, 16: 90354753.0, 17: 81195210.0, 18: 78077248.0, 19: 59128983.0, 2: 243199373.0,
			   20: 63025520.0, 21: 48129895.0, 22: 51304566.0, 3: 198022430.0, 4: 191154276.0, 5: 180915260.0,
			   6: 171115067.0, 7: 159138663.0, 8: 146364022.0, 9: 141213431.0, 23: 155270560.0, 24: 59373566.0}
from collections import defaultdict

###########################Rreading Gene#########################
###############################################################
genes = {}
genes = defaultdict(lambda:[],genes)
with open(args.gene , 'r') as f :
    for line in f :
        line = line.strip().split('\t')
        gene_name = line[-4]
        chromosome = line[2]
        strand = line[3]
        start = int(line[4])
        end = int(line[5])
        if chromosome.startswith('chr'):
            genes[chromosome].append([start, end , strand ,gene_name])
gene_interupt_neighbour = 20000

##########################################################################

a = {} # dict of chromosome number to dict of position to alignments

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
		a[chrom][j].append(i) #specify this alignment overlap which this cluster region
lines_to_write = []
with open(args.output, 'w')as file:
	file.write("#Type\tChromosome\tRefStartPos\tRefEndPos\tPredictedSize\tSupportsId\tSupportsCount\n")
	for k in a:# k = chromosome number
	    b = a[k] # b = ditc of cluster region to list of alignments
	    for kk in b: #kk = cluster region number
	        ans = {'insertion': [], 'deletion': []}
	        for i in b[kk]: #alignment overlap kk
	            alignment = i['Alignment']
	            q_id = i['QryContigID']
	            ref_id = i['RefContigID']
	            alignment = alignment.split(')')[:-1]
	            find_del = [] #
	            find_ins = []
	            for r in range(1, 4): #compare at most lengh window 4 started with 2 in alignment pair v
	                for v in range(0, len(alignment) - r):
	                    ref_1 = int(alignment[v].split(',')[0][1:]) # ref label number
	                    if kk * scale < ref[ref_id][ref_1] <= (kk + 1) * scale:
	                        q_1 = int(alignment[v].split(',')[1]) # query label number
	                        ref_2 = int(alignment[v + r].split(',')[0][1:])# ref2 label number
	                        q_2 = int(alignment[v + r].split(',')[1])# query2 label number
	                        ref_dist = abs(ref[ref_id][ref_1] - ref[ref_id][ref_2])
	                        q_dist = abs(mol[q_id][q_1] - mol[q_id][q_2])
	                        dif = ref_dist - q_dist
	                        if dif < -2000: #insertion with length greater than 2Kbp
	                        #check if this happen in not first 25 percent start and end of 
	                            if mol[q_id][min(q_1, q_2)] > min(0.25 * mol[q_id][len(mol[q_id])],10000) and mol[q_id][max(q_1, q_2)] < max(0.25 * mol[q_id][len(mol[q_id])],  mol[q_id][len(mol[q_id])]-10000):
	                                if not find_intersection(find_ins, (q_1, q_2)):
	                                    ans['insertion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
	                                    find_ins.append((q_1, q_2))

	                        if dif > 2000:
	                            if mol[q_id][min(q_1, q_2)] > min(0.25 * mol[q_id][len(mol[q_id])], 10000) and \
	                                    mol[q_id][
	                                        max(q_1, q_2)] < max(0.25 * mol[q_id][len(mol[q_id])],
	                                                             mol[q_id][len(mol[q_id])] - 10000):
	                                if not find_intersection(find_del, (q_1, q_2)):
	                                    ans['deletion'].append([q_id, ref_1, ref_2, q_1, q_2, dif])
	                                    find_del.append((q_1, q_2))

	        for i in ans:
	            if len(ans[i]) > 0: #chech is there any indel here in this region?
	                ans[i] = sorted(ans[i], key=lambda x: (x[1], x[2]))
	                selected = []#Selected SVid
	                i_j = (max(1, ans[i][0][1]), min(ans[i][0][2], len(ref[str(rev_contig[k])])))
	                s = 0 #Sum SV size
	                
	                for aln in ans[i]:
	                	#because it is sorted
	                    if aln[1] >= i_j[0] - 0 and aln[2] <= i_j[1] + 0:
	                        i_j = (min(i_j[0], aln[1]), max(i_j[1], aln[2]))
	                        selected.append(aln[0])
	                        s += aln[5]
	                    else:
	                        if len(selected) > 0:
	                            len_selected_before_set = len(selected)
	                            selected = list(set(selected))
	                            lines_to_write.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(i,k,ref[str(rev_contig[k])][i_j[0]],ref[str(rev_contig[k])][i_j[1]],float(float(s) / len_selected_before_set),','.join(selected),len(set(selected))))
	                            # file.write(i + '\t' + str(k) + '\t' + str(ref[str(rev_contig[k])][i_j[0]]) + '\t' + str(
	                            #     ref[str(rev_contig[k])][i_j[1]]) + '\t' + ','.join(selected) + '\t' + str(
	                            #     len(set(selected))) + '\t' + str(float(float(s) / len_selected_before_set)) + '\t' + str(
	                            #     float(float(s2) / len(selected))) + '\n')
	                        i_j = (max(1, aln[1] - 0), min(aln[2] + 0, len(ref[str(rev_contig[k])])))
	                        selected = [aln[0]]
	                        s = aln[5]
	                if len(selected) > 0:
	                    len_selected_before_set = len(selected)
	                    selected = list(set(selected))
	                    lines_to_write.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(i,k,ref[str(rev_contig[k])][i_j[0]],ref[str(rev_contig[k])][i_j[1]],float(float(s) / len_selected_before_set),','.join(selected),len(set(selected))))
	d = {} #dict molecule ID to partial alignment of this molecule
	#part2 looking for indel in 2 separate partial alignments
	d = defaultdict(lambda: [], d)
	with open(args.alignment, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				line2 = line.strip().split('\t')
				id = int(line2[1])
				d[id].append(line2)
	for k in d.keys(): #k = molecule ID
		for i in range(len(d[k])):
			q = d[k][i] #partial alignment 1
			align = q[-1]
			chrom = q[2]
			s = int(align.split(',')[1].split(')')[0]) #molecule start label number
			s_r = int(align.split(',')[0][1:]) #ref start label number
			s_ref = q[5] #reference start coordiante
			e = int(
				align.split(',')[-1].split(')')[0])  #molecule end label number
			e_r =int(align.split(',')[-2].split('(')[1]) #ref end label number
			e_ref = q[6] #reference end coordiante
			direction = q[7]
			if direction == '+':
				next_s = s - 1 # next molecule label expected in start
				next_e = e + 1 # next molecule label expected in end
				next_s_q = mol[str(k)][max(1, next_s)] # next molecule coordinate expected in start
				next_e_q = mol[str(k)][min(next_e, len(mol[str(k)]))] # next molecule coordinate expected in end
			else:
				next_s = s + 1
				next_e = e - 1
				next_s_q = mol[str(k)][min(next_s, len(mol[str(k)]))]
				next_e_q = mol[str(k)][max(1, next_e)]
			for j in range(len(d[k])):
				if j != i:
					x1 = 0
					x2 = 0
					p = d[k][j] #another partial alignment
					if p[7] == direction and p[2] == chrom:
						#check have overlap with expected one
						if min(float(p[3]),float(p[4])) <= next_s_q <= max(float(p[3]),float(p[4])) or min(float(p[4]),float(p[3])) <= next_e_q <= max(float(p[3]),float(p[4])):
							align_p = p[-1].split('(')[1:]
							selected_pair = align_p[:3]+align_p[-3:]
							for pair in selected_pair:
								if str(next_s) == pair.split(',')[1][:-1]:
									s_ref_p = ref[str(chrom)][int(pair.split(',')[0])]
									#check if distance is greater than 2KBp and not to faar
									if abs(float(s_ref) - float(s_ref_p)) - abs(next_s_q - mol[str(k)][s]) >2000 and  abs(float(s_ref) - float(s_ref_p)) < 300000:
										lines_to_write.append("deletion\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(contig_id_map[int(chrom)],s_ref_p,s_ref,abs(float(s_ref) - float(s_ref_p)) - abs(next_s_q - mol[str(k)][s]),k,1))
									x1 = 1
								if str(next_e) == pair.split(',')[1][:-1]:
									e_ref_p = ref[str(chrom)][int(pair.split(',')[0])]
									if abs(float(e_ref) - float(e_ref_p)) - abs(next_e_q - mol[str(k)][e]) >2000  and abs(float(e_ref) - float(e_ref_p)) < 300000: 
										lines_to_write.append("deletion\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(contig_id_map[int(chrom)],e_ref,e_ref_p,abs(float(e_ref) - float(e_ref_p)) - abs(next_e_q - mol[str(k)][e]),k,1))
									x2 = 1
								if x1 + x2 == 2:
									break
	lines_to_write = set(lines_to_write)
	# for line_t in lines_to_write:
	# 	file.write(line_t)
	lines_to_write = sorted(lines_to_write, key=lambda x: (x.strip().split('\t')[0], int(x.strip().split('\t')[1]), float(x.strip().split('\t')[2]), float(x.strip().split('\t')[3])))
	d = {}
	d = defaultdict(lambda: defaultdict(list), d)
	for line_t in lines_to_write:
		c = int(line.split('\t')[1])
		sv_type = line.split('\t')[0]
		d[c][sv_type].append(line_t)
	for k in d:
		for sv_type in d[k]:
			a  = d[k][sv_type]
			if len(a)> 0:
				w = []
				for i in range(len(a)-1):
					aln = a[i]
					nex_aln = a[i+1]
					aln_2 = aln.strip().split('\t')
					nex_aln_2 = nex_aln.strip().split('\t')
					if abs(float(aln_2[3]) - float(nex_aln_2[2])) < 30000:
						nex_aln_2[2] = aln_2[2]
						nex_aln_2[3] = str(max(float(nex_aln_2[3]),float(aln_2[3])))
						nex_aln_2[5] +=','+aln_2[5]
						nex_aln_2[5] = ','.join(list(set(nex_aln_2[5].split(','))))
						nex_aln_2[6] = str(len(nex_aln_2[5].split(',')))
						nex_aln_2[4] = str(float(nex_aln_2[4]) + float(aln_2[4]))
						a[i+1] = '\t'.join(nex_aln_2)+ '\n'
					elif (abs(float(aln_2[2]) - float(nex_aln_2[2])) < 30000 and abs(float(aln_2[3]) - float(nex_aln_2[3])) < 30000 ):
						nex_aln_2[2] = aln_2[2]
						nex_aln_2[3] = str(max(float(nex_aln_2[3]),float(aln_2[3])))
						nex_aln_2[5] +=','+aln_2[5]
						nex_aln_2[5] = ','.join(list(set(nex_aln_2[5].split(','))))
						nex_aln_2[6] = str(len(nex_aln_2[5].split(',')))
						nex_aln_2[4] = str(max(float(nex_aln_2[4]) , float(aln_2[4])))
						a[i+1] = '\t'.join(nex_aln_2)+ '\n'
					else:
						w.append(a[i])
				w.append(a[-1])
				for i in w:
					i = i.strip().split('\t')
					chrom = i[1]
					pos_start = float(i[2])
					pos_end = float(i[3])
					if str(chrom) =='23':
						chrom = 'X'
					if str(chrom) =='24':
						chrom = 'Y'
					chrom = 'chr'+ str(chrom)
					gene_interupt = []
					genes_chr = genes[chrom]
					for g in genes_chr:
						if max(pos_start, g[0]) < min (pos_end,g[1]):
							gene_interupt.append(g[-1])
					gene_interupt = list(set(gene_interupt))
					i.append(','.join(gene_interupt))
					file.write('\t'.join(i)+'\n')
