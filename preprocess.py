import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="query directory", required=True)
parser.add_argument("-o", "--output", help="Output dir for processed contigs", required=True)
parser.add_argument("-m", "--maxsize", help="Maximum size for each molecule. Default is 100", required=False)
args = parser.parse_args()
max_len = 100
if str(args.maxsize)!='None':
	max_len = int(args.maxsize)

long_lines = []
max_id = 0
map_new_id = {}
with open(args.query, 'r') as f:
    with open(args.output+'.cmap', 'w') as g:
        for line in f:
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                if int(line2[0]) > max_id:
                    max_id =int(line2[0])
                if int(line2[2]) < max_len+1:
                    g.write(line)
                else:
                    long_lines.append(line)
        id = 0
        number_of_split = 1
        last_mol_len = 0
        mol_len = 0
        all_len = 0
        for i,line in enumerate(long_lines):
            line2 = line.strip().split('\t')
            if int(line2[0]) !=id:
                id = int(line2[0])
                max_id = max_id + 1
                number_of_split = int(line2[2]) // max_len
                map_new_id[max_id] = [id,'0','0']
                counter_split = 1
                mol_len = long_lines[i+max_len-1].split('\t')[5]
                all_len = int(line2[2])
                last_mol_len = 0
            line2[0] = str(max_id)
            if counter_split > number_of_split:
                line2[2] = str(int(line2[2]) % max_len)
            else:
                line2[2]= str(max_len)
            line2[5] = str(round(float(line2[5]) - last_mol_len,1))
            line2[3] = str(int(line2[3]) % max_len)
            if line2[3] == '0':
                line2[3] = str(max_len)
            line2[1] = mol_len
            if int(line2[3]) % max_len == 0:
                max_id +=1
		last_mol_len =last_mol_len +  float(mol_len)
                map_new_id[max_id] = [id, str(counter_split * max_len),str(round(last_mol_len,1))]
                counter_split += 1
                if counter_split > number_of_split:
                    mol_len = str(round(float(long_lines[i + all_len %max_len].split('\t')[5]) - float(last_mol_len),1))
                else:
                    mol_len = str(round(float(long_lines[i+max_len].split('\t')[5]) - float(last_mol_len),1))
                if not(all_len % max_len ==0 and counter_split > number_of_split):
                    g.write('\t'.join(line2))
                    g.write('\n')        
                    line2[4] = '0'
                    line2[3] = str(int(line2[3])+1)

            if line2[4] == '0':
                line2[5] = line2[1]
            g.write('\t'.join(line2))
            g.write('\n')


with open(args.output+'_dic','w') as f:
    for k in map_new_id:
        f.write(str(k) + '\t'+ str(map_new_id[k][0])+'\t'+map_new_id[k][1]+'\t' +map_new_id[k][2]+'\n')

