#high_conf = '/home/sraeisid/NA12878_svs.vcf'
high_conf = '/home/sraeisid/high_conf'
from collections import defaultdict

SV={}
for k in range(1,25):
    SV[k] = []
    # SV[k] = defaultdict(lambda: [], SV[k])
con = 0
with open (high_conf, 'r')as f:
	for line in f :
		if not line.startswith('#'):
			line2 = line.split('\t')
			chrom = line2[0]
			if chrom == 'X':
				chrom = 23
			elif chrom == 'Y':
				chrom = 24
			elif chrom == 'MT':
				continue
			else:
				chrom = int(chrom)
			start_pos = int(line2[1])
			type_sv = line2[4][1:-1]
			end_pos = int(line2[7].split(';')[0].split('=')[1])
			start_pos,end_pos = min(start_pos,end_pos),max(start_pos,end_pos)
			if 30000>=abs(end_pos - start_pos) > 2000 and type_sv!='INV':
				con +=1
				SV[chrom].append([type_sv,start_pos,end_pos])


NA_dir = '/nucleus/projects/sraeisid/OMFilter_result/NA12878/test_assemble/del2_cluster'
d_s = []
d_d = []
all = 0
correct = 0
del_c = 0
in_c = 0
with open(NA_dir , 'r') as f:
	for line in f:
		line2 = line.split('\t')
		sv_type= line2[0]
		chrom = int(line2[1])
		start_pos, end_pos = min(float(line2[2]),float(line2[3])), max(float(line2[2]),float(line2[3]))
		all +=1
		find = 0
		for a  in SV[chrom]:
			if sv_type =='deletion':
				 if a[0]=='DEL' and a[1]<=end_pos and start_pos <= a[2]:
#				 	print(a)
#					print(line.strip())
					SV[chrom].remove(a)
					d_s.append(float(line2[7].strip()))
				 	correct+=1
				 	del_c +=1
				 	find = 1
				 	break
			else:
				if a[0]=='INS' and a[1]<=end_pos and start_pos <= a[2]:
					SV[chrom].remove(a)
					d_s.append(float(line2[7].strip()))
#					print(a)
					correct+=1
					in_c +=1
					break
		d_d.append(float(line2[7].strip()))
#		if find ==0:
 #               	print(line.strip())
print('Number in Fandom',all)
print('Number find in gold standard',correct)

print(float(float(correct)/all))
print('Number detectable in gold standard',con)
print('find/detectable',float(float(correct)/con))
print(sum(d_s)/len(d_s))
print(sum(d_d)/len(d_d))
print('Insertion', in_c)
print('deletion', del_c)
# for chrom in SV.keys():
# 	for j in SV[chrom]:
# 		print(int(j[2]-j[1]))






