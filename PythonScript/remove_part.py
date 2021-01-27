
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
                d.append(line[1])
    d = list(set(d))
    return d
def remove_partial(d):
    with open(args.full,'r') as f:
        with open(args.output+'.xmap','w')as g :
            for line in f :
                if line.startswith('#'):
                    g.write(line)
                else:
                    line2 = line.strip().split('\t')
                    if not line2[1] in d:
                        g.write(line)
                    
if __name__ == '__main__':
    d = parse_partial_alignments()
    remove_partial(d)