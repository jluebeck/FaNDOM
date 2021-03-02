from collections import defaultdict
from scipy.stats import entropy
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output Xmap file", required=True)
parser.add_argument("-r", "--ref", help="Directory to Reference", required=False)
args = parser.parse_args()
d = {}
d = defaultdict(lambda: [], d)
perv = 0
with open(args.ref, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            id = int(line[0])
            pos = float(line[5])
            d[id].append(pos)
low_complexity ={}
low_complexity = defaultdict(lambda: [], low_complexity)
for k in d.keys():
    distances = d[k]
    for i in range(0, len(distances) - 6):
        dis = distances[i:i + 6]
        b = []
        for j in range(0, len(dis) - 1):
            b.append(dis[j + 1] - dis[j])
        a = entropy(b, base = 5)
        if abs(a - 1) < 0.005:
            for j in range(0, len(dis) - 1):
                low_complexity[k].append(i+j+1)
            low_complexity[k].append(i + 6)
for k in low_complexity.keys():
    low_complexity[k] = sorted(list(set(low_complexity[k])))
d = {}
kk = 1;
from collections import defaultdict
d = defaultdict(lambda:[],d)
with open(args.output+'.xmap', 'w') as g:
    with open(args.input, 'r') as f:
        for line in f:
            if not line.startswith('#'):
		# print(kk)
		kk+=1
                line2 = line.strip().split('\t')
                id = int(line2[1])
                conf = float(line2[8])
                a1 = float(line2[3])
                a2 = float(line2[4])
                align = line2[13]
		a = align.split((')('))
		a[0] = a[0][1:]
		count = 0
		for i in a:
			r = i.split(',')[0]
			r = int(r)
			if r in  low_complexity[int(line2[2])]:
				count+=1
			else:
				count = count - 0.2 
		if count > 5:
			continue
                b1 = int(align.split(')')[0].split(',')[1])
                b2 = int(align.split(')')[-2].split(',')[1])
                align = align.split(',')
                separate_lines = len(align) - 1
                score = conf / separate_lines
                if abs(a1 - a2) > 30000 and abs(b1 - b2) > 8 and score > 4000:
			#if not (score < 5000 and line2[2]=='25'):
			d[id].append((score, line))
            else:
                g.write(line)

    for k in d.keys():
	# print('ali')
        d[k] = sorted(d[k], reverse=True)[:300]

    for k in d.keys():
	# print('reaz')
	# print(k)
        deleted = []
        for i in range(len(d[k])):
            a = d[k][i]
            a1 = a[1].split('\t')
            for j in range(i + 1, len(d[k])):
                b = d[k][j]
                b1 = b[1].split('\t')
                if min(float(a1[3]), float(a1[4])) <= min(float(b1[3]), float(b1[4])) and max(float(a1[3]),
                                                                                              float(a1[4])) >= max(
                    float(b1[3]), float(b1[4])):
                    deleted.append(b)
        d[k] = [x for x in d[k] if not x in deleted]
    for k in d.keys():
        for a in d[k]:
            g.write(a[1])

