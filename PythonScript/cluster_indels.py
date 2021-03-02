import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Indels input", required=True)
parser.add_argument("-o", "--output", help="Clustered Indels Output", required=True)
args = parser.parse_args()


def parse_inde():
    d = {}
    d = defaultdict(lambda: [], d)
    with open(args.input, 'r') as f:
        for line in f:
            chromosome = int(line.split('\t')[1])
            d[chromosome].append(line)
    return d




if __name__ == '__main__':
    d = parse_inde()
    with open(args.output, 'w') as f:
        for chromosome in d:
            a = d[chromosome]
            if len(a) > 0:
                w = []
                for i in range(len(a) - 1):
                    aln = a[i]
                    nex_aln = a[i + 1]
                    aln_2 = aln.strip().split('\t')
                    nex_aln_2 = nex_aln.strip().split('\t')

                    if abs(float(aln_2[3]) - float(nex_aln_2[2])) < 30000:
                        nex_aln_2[2] = aln_2[2]
                        nex_aln_2[3] = str(max(float(nex_aln_2[3]), float(aln_2[3])))
                        nex_aln_2[5] += ',' + aln_2[5]
                        nex_aln_2[4] = str(float(nex_aln_2[4]) + float(aln_2[4]))
                        a[i + 1] = '\t'.join(nex_aln_2) + '\n'
                    elif (abs(float(aln_2[2]) - float(nex_aln_2[2])) < 30000 and abs(
                            float(aln_2[3]) - float(nex_aln_2[3])) < 30000):
                        nex_aln_2[2] = aln_2[2]
                        nex_aln_2[3] = str(max(float(nex_aln_2[3]), float(aln_2[3])))
                        nex_aln_2[5] += ',' + aln_2[5]
                        nex_aln_2[4] = str(max(float(nex_aln_2[4]), float(aln_2[4])))
                        a[i + 1] = '\t'.join(nex_aln_2) + '\n'
                    else:
                        w.append(a[i])
                w.append(a[-1])
                for i in w:
                    f.write(i)