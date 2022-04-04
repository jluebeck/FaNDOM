from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output Xmap file", required=True)
#parser.add_argument("-r", "--ref", help="Directory to Reference", required=False)
args = parser.parse_args()

def nov_overlap(a1,a2,a3,a4):
    if max(a1,a3)== a1:
        return a2 - a4
    else:
        return a3 - a1

def span_alignment(a1,a2,a3,a4):
    minimum_nonoverlap_length = 40000
    a1,a2 = min(a1,a2), max(a1,a2)
    a3,a4 = min(a3,a4), max(a3,a4)
    span_area = min(a2,a4) - max(a1,a3) 
    answer_div = float(span_area/(a4-a3))
    if answer_div >= 0.8:
        return answer_div
    elif answer_div >= 0.5:
        if nov_overlap(a1,a2,a3,a4) >= minimum_nonoverlap_length:
            return 1
    return 0 
d = {}
# kk = 1;
d = defaultdict(lambda: [], d)
with open(args.output + '.xmap', 'w') as g:
    with open(args.input, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line2 = line.strip().split('\t')
                id = int(line2[1])
                conf = float(line2[8])
                a1 = float(line2[3])
                a2 = float(line2[4])
                align = line2[13]
                align = align.split(',')
                separate_lines = len(align) - 1
                score = conf / separate_lines
                if abs(a1 - a2) > 40000 and ( (12> abs(separate_lines) > 9 and score > 5000) or (16> abs(separate_lines) > 11 and score > 4500 ) or ( abs(separate_lines) >15 and score >4000 ) ):
                    # if not (score < 5000 and line2[2]=='25'):
                    d[id].append((score, line))
            else:
                g.write(line)

    for k in d.keys():
        d[k] = sorted(d[k], reverse=True)[:300]

    for k in d.keys():

        # print(k)
        deleted = []
        for i in range(len(d[k])):
            a = d[k][i]
            a1 = a[1].split('\t')
            for j in range(i + 1, len(d[k])):
                b = d[k][j]
                b1 = b[1].split('\t')
                if span_alignment(float(a1[3]), float(a1[4]), float(b1[3]), float(b1[4])) >=0.8:
                    deleted.append(b)
                elif min(float(b1[3]), float(b1[4])) - 2000 <= min(float(a1[3]), float(a1[4])) and max(float(b1[3]),
                                                                                                float(b1[4])) + 2000 >= max(
                    float(a1[3]), float(a1[4])):
                                                                                                if abs(float(b1[4]) -float(b1[3])) - abs(float(a1[4]) - float(a1[3])) > 2000 :
                                                                                                    deleted.append(a)

        d[k] = [x for x in d[k] if not x in deleted]
    for k in d.keys():
        for a in d[k]:
            g.write(a[1])
