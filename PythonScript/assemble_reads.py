def join_string(main , next , direction):
    if next in main:
        return main
    if direction == '+':
        while True:
            a = next.split(')')[0]
            a = a + ')'
            if a in main:
                break
            else:
                next = next.split(')')[1:]
                next = ')'.join(next)
        while True:
            a = main.split(')')[-2]
            a = a+ ')'
            if a in next:
                break
            else:
                main =main[:main.rfind('(')]
        if next in main:
            return  main
        main = main.split(')')[:-4]
        main = ')'.join(main)+ ')'
        find =1
        i = 1
        while find ==1 :
            if next[:i] in main:
                i+=1
            else:
                find = 0
        return main + next[i-1:]
    else:
        while True:
            a = next.split(')')[-2]
            a = a+ ')'
            if a in main:
                break
            else:
                next =next[:next.rfind('(')]
        while True:
            a = main.split(')')[0]
            a = a + ')'
            if a in next:
                break
            else:
                main = main.split(')')[1:]
                main = ')'.join(main)
        if next in main:
            return  main
        main = main.split(')')[3:]
        main = ')'.join(main)
        find = 1
        i = 1
        while find ==1:
            if next[len(next)-i:] in main:
                i = i+ 1
            else:
                find = 0
        return next[:len(next) -i+1] + main
from collections import defaultdict
d = {}
d = defaultdict(lambda:[],d)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for Xmap file", required=True)
args = parser.parse_args()

with open(args.output+'_assembled.xmap','w') as g:
    with open(args.input,'r')as f:
        for line in f:
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                d[line2[1]].append(line2)
    for k in d :
	print(k)
        d[k] = sorted(d[k], key=lambda x: min(float(x[3]),float(x[4])))
        i = 1
        main_aln = d[k][0]
        while i < len(d[k]):
            find = 0
            nex_aln = d[k][i]
            if main_aln[2] == nex_aln[2] and main_aln[7] == nex_aln[7]:
                qq = max(float(main_aln[3]), float(main_aln[4]))
                ww = min(float(nex_aln[3]), float(nex_aln[4]))
                if max(float(main_aln[3]), float(main_aln[4])) > min(float(nex_aln[3]), float(nex_aln[4])) and len(set(main_aln[-1].split(')')) & set(nex_aln[-1].split(')'))) > 1:
                    if max(float(main_aln[5]), float(main_aln[6])) > min(float(nex_aln[5]), float(nex_aln[6])) and min(float(main_aln[5]), float(main_aln[6])) < min(float(nex_aln[5]), float(nex_aln[6])) and main_aln[7] == '+':
                        if main_aln[7] == '+' :
                            main_aln[4] = nex_aln[4]
                            main_aln[6] = nex_aln[6]
                            main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '+')
                        else:
                            main_aln[3] = nex_aln[3]
                            main_aln[5] = nex_aln[5]
                            main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '-')
                        i+=1
                        find = 1
                    if max(float(nex_aln[5]), float(nex_aln[6])) > min(float(main_aln[5]), float(main_aln[6])) and min(float(nex_aln[5]), float(nex_aln[6])) < min(float(main_aln[5]), float(main_aln[6])) and main_aln[7] == '-':
                        if main_aln[7] == '+' :
                            main_aln[4] = nex_aln[4]
                            main_aln[6] = nex_aln[6]
                            main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '+')
                        else:
                            main_aln[3] = nex_aln[3]
                            main_aln[5] = nex_aln[5]
                            main_aln[-1] = join_string(main_aln[-1], nex_aln[-1], '-')
                        i+=1
                        find = 1

            if find == 0:
                g.write('\t'.join(main_aln))
                g.write('\n')
                main_aln = nex_aln
                i+=1
        g.write('\t'.join(main_aln))
        g.write('\n')




