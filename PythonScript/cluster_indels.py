import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Indels input", required=True)
parser.add_argument("-o", "--output", help="Clustered Indels Output", required=True)
args = parser.parse_args()

d = {}
for i in range(25):
	d[i] = []
with open(args.input,'r') as f:
	for line in f :
		c = int(line.split('\t')[1])
		d[c].append(line)
with open(args.output,'w') as f:
	for k in d:
		a  = d[k]
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
					nex_aln_2[4] +=','+aln_2[4]
					nex_aln_2[6] = str(float(nex_aln_2[6]) + float(aln_2[6]))
					a[i+1] = '\t'.join(nex_aln_2)+ '\n'
				elif (abs(float(aln_2[2]) - float(nex_aln_2[2])) < 30000 and abs(float(aln_2[3]) - float(nex_aln_2[3])) < 30000 ):
					nex_aln_2[2] = aln_2[2]
					nex_aln_2[3] = str(max(float(nex_aln_2[3]),float(aln_2[3])))
					nex_aln_2[4] +=','+aln_2[4]
					nex_aln_2[6] = str(max(float(nex_aln_2[6]) , float(aln_2[6])))
					a[i+1] = '\t'.join(nex_aln_2)+ '\n'
				else:
					w.append(a[i])
			w.append(a[-1])
			for i in w:
				f.write(i)




