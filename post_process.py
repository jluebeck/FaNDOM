import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fandom", help="Path to FaNDOM alignment", required=True)
parser.add_argument("-d", "--dic", help="Path to Dictionary made by preprocess ", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
args = parser.parse_args()

d = {}
with open(args.dic,'r') as f:
    for line in f :
        line = line.strip().split('\t')
        d[int(line[0])] = [line[1], int(line[2]), float(line[3])]

with open(args.fandom,'r') as f :
    with open(args.output+'_post_process.xmap', 'w') as g:
        for line in f :
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                if int(line2[0]) in d.keys():
                    c = d[int(line2[0])][1]
		    shift = d[int(line2[0])][2]
                    line2[0] = d[int(line2[0])][0]
                    alignment = line2[-1]
                    alignment = alignment.split(')')[:-1]
                    new_align = ''
                    for i in alignment:
                        b = i.split(',')
                        new_align =new_align +  b[0]+','+str(int(b[1])+c)+')'
                    line2[-1] = new_align
                    line2[3] = str(round(float(line2[3]) + shift,1))
                    line2[4] = str(round(float(line2[4]) + shift,1))
                g.write('\t'.join(line2))
                g.write('\n')

