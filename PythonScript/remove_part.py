from bisect import bisect_left  
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--partial", help="Partial Alignment", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-f", "--full", help="Full Alignment", required=True)

args = parser.parse_args()
def parse_partial_alignments():
    d = []
    with open(args.partial,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split('\t')
                d.append(int(line[1]))
    d = list(set(d))
    d.sort()
    return d
def remove_partial(d):
    with open(args.full,'r') as f:
        with open(args.output+'.xmap','w')as g :
            for line in f :
                if line.startswith('#'):
                    g.write(line)
                else:
                    line2 = line.strip().split('\t')
                    id = int(line2[1])
                    index = bisect_left(d, id)
                    if index == len(d):
                        g.write(line)
                    elif d[index]!=id:
                        g.write(line)
                    
if __name__ == '__main__':
    d = parse_partial_alignments()
    remove_partial(d)