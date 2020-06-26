def parse_cmap(cmapf, keep_length=False):
    cmaps = {}
    # contigCovs = {}
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


def parse_bnx(bnxF, keep_length=False):
    moleculeD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                elif line.startswith('1'):
                    if keep_length:
                        moleculeD[currKey] = [float(x) for x in fields[1:]]
                    else:
                        moleculeD[currKey] = [float(x) for x in fields[1:-1]]

    return moleculeD


import argparse
import random

sampleSize = 500

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="query directory", required=True)
parser.add_argument("-r", "--ref", help="reference directory", required=True)
parser.add_argument("-f", "--fandom", help="FaNDOM directory", required=True)
parser.add_argument("-o", "--output", help="Output directory", required=True)
parser.add_argument("-t", "--thread", help="Output directory", required=True)
args = parser.parse_args()

if args.query.endswith('.bnx'):
    query = parse_bnx(args.query,True)
elif args.query.endswith('.cmap'):
    query = parse_cmap(args.query,True)
else:
    print('Query file is not in cmap or bnx format')
    return
if args.fandom.endswith('/'):
    fandom = args.fandom
else:
    fandom = args.fandom+'/'

if !args.ref.endswith('.cmap'):
    print('Reference should be in cmap format')
    return
if args.output.endswith('/'):
    output = args.output[:-1]
else:
    output = args.output
selected = [list(query.keys())[i] for i in random.sample(range(0, len(query)), sampleSize)]
import os
scaledSample = {}
result = {}
for i in range(0, 11):
    with open(output+'/a.cmap','w') as f :
        for k in selected:
            scaledSample[k] = query[k].copy()
            if args.query.endswith('.bnx'):
                scaledSample[k] = {i+1:j for i,j in enumerate(scaledSample[k])}
            for index in scaledSample[k].keys():
                scaledSample[k][index] =float(scaledSample[k][index] * (100-i))/100
            length = scaledSample[k][list(scaledSample[k])[-1]]
            number_label = index
            for index in scaledSample[k]:
                if index ==number_label:
                    f.write(k + '\t' + str(length) + '\t' + str(number_label) + '\t' + str(index) + '\t0\t' + str(
                        scaledSample[k][index]) + '\t0.0\t1.0\t1.0\t17.0000\t0.0000\n')
                else:
                    f.write(k + '\t'+str(length)+'\t'+str(number_label)+'\t'+str(index)+'\t1\t'+str(scaledSample[k][index])+'\t0.0\t1.0\t1.0\t17.0000\t0.0000\n')

    os.system(fandom+'./FaNDOM'+' -t='+args.thread+' -r='+args.ref+' -q='+output+'/a.cmap -sname='+output+'/test '+ '-outfmt=xmap')
    sum_score = 0
    with open(output+'/test.xmap','r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split('\t')
                sum_score += float(line[8])
    result[i] = sum_score
print(result)
os.system('rm '+output+'/a.cmap')
max_index = int(max(result,key=result.get))

print('Max score with rescale factor: ', max_index)
if args.query.endswith('.bnx'):
    with open(output+'/rescaled.bnx','w') as f:
        with open(args.query,'r') as g:
            for line in g:
                if line.startswith('0\t'):
                    line = line.strip().split('\t')
                    line[2] = str((float(line[2]) * (100-max_index))/100)
                    line = '\t'.join(line) + '\n'
                elif line.startswith('1\t'):
                    line = line.strip().split('\t')
                    for i in range(1,len(line)):
                        line[i] = str((float(line[i]) * (100-max_index))/100)
                    line = '\t'.join(line) + '\n'
                f.write(line)

if args.query.endswith('.cmap'):
    with open(output + '/rescaled.cmap', 'w') as f:
        with open(args.query, 'r') as g:
            for line in g:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    line[1] = str((float(line[1]) * (100-max_index))/100)
                    line[5] = str((float(line[5]) * (100 - max_index)) / 100)
                    line = '\t'.join(line) + '\n'
                f.write(line)
