import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Xmap file", required=True)
parser.add_argument("-o", "--output", help="Output Xmap file", required=True)

args = parser.parse_args()

contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}

if __name__ == '__main__':
	with open(args.output+'.xmap' ,'w') as f:
		with open(args.input,'r') as g:
			for line in g:
				if not line.startswith('#'):
					line = line.strip().split('\t')
					line[2] = str(contig_id_map[int(line[2])])
					line = '\t'.join(line)
					f.write(line)
					f.write('\n')
				else:
					f.write(line)

