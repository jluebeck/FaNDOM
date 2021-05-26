import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input reference genome in cmap format", required=True)
parser.add_argument("-o", "--output", help="Output directory in cmap", required=True)
args = parser.parse_args()
def merge_distance(reference, tolerance, output):
    d ={}
    all = []
    with open(reference) as infile:
        with open(output, 'w') as the_file:
            first = 1
            save_line = ''
            index = 1
            for line in infile:
                if line.startswith("#"):
                        the_file.write(line)
                elif not line.startswith("#"):
                    fields = line.rstrip().rsplit()
                    if first:
                        first = 0
                        save_line = fields
                    else:
                        if int(save_line[0]) == int(fields[0]):
                            if abs(float(save_line[5]) - float(fields[5])) > tolerance:
                                save_line[3] = index
                                all.append(save_line)
                                index += 1
                                save_line = fields
                            else:
                                save_line[5] = (float(save_line[5]) + float(fields[5])) / 2
                        else:
                            
                            save_line[3] = index
                            d[save_line[0]] = int(save_line[3])-1
                            all.append(save_line)
                            index = 1
                            save_line = fields
            save_line[3] = index
            d[save_line[0]] = int(save_line[3])-1
            all.append(save_line)
            for s in all:
                s[2] = str(d[s[0]])
                the_file.write('\t'.join([str(x) for x in s]))
                the_file.write('\n')


                
ref = args.input
T = 800
output = args.output

merge_distance(ref,T,output)
