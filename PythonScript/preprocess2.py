import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="query directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-m", "--maxsize", help="Maximum size for each molecule. Default is 25", required=False)
args = parser.parse_args()
max_len = 25
if str(args.maxsize)!='None':
	max_len = int(args.maxsize)
const = 50
long_lines = []
max_id = 0
map_new_id = {}
l = []

def sep(lines, last_mol_len, counter_split, last):
    global map_new_id
    global max_id
    global l
    id = int(lines[0].split('\t')[0])
    # with open('/Users/siavashraeisidehkordi/Desktop/University/Vineet_Research/indel/out2.cmap', 'a') as g:
    last_mol_len = float(last_mol_len)
    counter_split = counter_split - 1
    map_new_id[max_id] = [id, str(counter_split * max_len), str(round(last_mol_len, 1))]
    mol_len = str(round(float(lines[-1].split('\t')[5]) - float(last_mol_len), 1))
    if last == 0:
        for i, line in enumerate(lines):
            line2 = line.strip().split('\t')
            line2[0] = str(max_id)
            line2[2] = str(len(lines))
            line2[5] = str(round(float(line2[5]) - last_mol_len, 1))
            line2[3] = str(i + 1)
            line2[1] = mol_len
            # a = a + '\t'.join(line2)
            # a = a + '\n'
            # g.write('\t'.join(line2))
            l.append('\t'.join(line2))
            # g.write('\n')
        line2[4] = '0'
        line2[3] = str(int(line2[3]) + 1)
        line2[5] = line2[1]
        # a = a + '\t'.join(line2)+'\n'
        # g.write('\t'.join(line2))
        # g.write('\n')
        l.append('\t'.join(line2))
    else:
        for i, line in enumerate(lines):
            line2 = line.strip().split('\t')
            line2[0] = str(max_id)
            line2[2] = str(len(lines) - 1)
            line2[5] = str(round(float(line2[5]) - last_mol_len, 1))
            line2[3] = str(i + 1)
            line2[1] = mol_len
            # a = a + '\t'.join(line2)
            # a = a + '\n'
            # g.write('\t'.join(line2))
            # g.write('\n')
            l.append('\t'.join(line2))
    max_id += 1
    # g.close()


def seprate_mol(lines, const, max_line):
    number_of_split = len(lines) // max_line
    i = 1
    if number_of_split == 1:
        if len(lines) % max_line > 40:
            sep(lines[:max_line + const], 0, 1, 0)
            last_mol_len = lines[max_line - 1].split('\t')[5]
            sep(lines[max_line:], last_mol_len, 2, 1)
        else:
            sep(lines, 0, 1, 1)
    else:
        if len(lines) % max_line > 40:
            while i <= number_of_split:
                if i == 1:
                    last_mol_len = 0
                else:
                    last_mol_len = lines[(i - 1) * max_line - 1].split('\t')[5]
                sep(lines[(i - 1) * max_line:i * max_line + const], last_mol_len, i, 0)
                i += 1
            last_mol_len = lines[(i - 1) * max_line - 1].split('\t')[5]
            sep(lines[number_of_split * max_line:], last_mol_len, i, 1)
        else:
            while i < number_of_split:
                if i == 1:
                    last_mol_len = 0
                else:
                    last_mol_len = lines[(i - 1) * max_line - 1].split('\t')[5]
                sep(lines[(i - 1) * max_line:i * max_line + const], last_mol_len, i, 0)
                i += 1
            if i == 1:
                last_mol_len = 0
            else:
                last_mol_len = lines[(i - 1) * max_line - 1].split('\t')[5]
            sep(lines[(i - 1) * max_line:], last_mol_len, i, 1)


with open(args.query,'r') as f:
    with open(args.output+'.cmap', 'w') as g:
        for line in f:
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                if int(line2[0]) > max_id:
                    max_id = int(line2[0])
                if int(line2[2]) < max_len + 1:
                    g.write(line)
                else:
                    long_lines.append(line)
        id = long_lines[0].split('\t')[0]
        max_id += 1
        lines = []
        for i in long_lines:
            if id == i.split('\t')[0]:
                lines.append(i)
            else:
                id = i.split('\t')[0]
                seprate_mol(lines, const, max_len)
                lines = [i]
        seprate_mol(lines, const, max_len)

with open(args.output+'_dic', 'w') as f:
    for k in map_new_id:
        f.write(str(k) + '\t' + str(map_new_id[k][0]) + '\t' + map_new_id[k][1] + '\t' + map_new_id[k][2] + '\n')
with open(args.output+'.cmap', 'a') as g:
    for i in l:
        g.write(i)
        g.write('\n')
