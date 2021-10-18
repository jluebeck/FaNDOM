import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input_File", required=True)
parser.add_argument("-o", "--output", help="Output_dir", required=True)
parser.add_argument("-n", "--name", help="Input Name", required=True)
args = parser.parse_args()

input_dir = args.input
output_dir = args.output
name = args.name
header = ''
print_line = ''
i = 0
counter= 1
if input_dir.endswith('bnx'):
    with open(input_dir, 'r') as f:
        for line  in f :
            if line.startswith('#'):
                header = header + line
            else:
                i=i+1
                print_line = print_line + line
                if i == 500000:
                    with open(output_dir+'/'+name+'_'+str(counter)+'.bnx', 'w')as g:
                        g.write(header)
                        g.write(print_line)
                        i=0
                        counter+=1
                        print_line = ''
        # print(header)
        with open(output_dir + '/' + name + '_' + str(counter) + '.bnx', 'w')as g:
            g.write(header)
            g.write(print_line)

elif input_dir.endswith('cmap'):
    mold_id = '0'
    with open(input_dir, 'r') as f:
        for line  in f :
            if line.startswith('#'):
                header = header + line
            else:
                if line.split('\t')[0]!=mold_id:
                    mold_id = line.split('\t')[0]
                    i=i+1
                if i == 125000:
                    with open(output_dir+'/'+name+'_'+str(counter)+'.cmap', 'w')as g:
                        g.write(header)
                        g.write(print_line)
                        i=0
                        counter+=1
                        print_line = ''

                print_line = print_line + line
                
        # print(header)
        with open(output_dir + '/' + name + '_' + str(counter) + '.cmap', 'w')as g:
            g.write(header)
            g.write(print_line)    



