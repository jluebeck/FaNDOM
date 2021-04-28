input_dir = '/nucleus/projects/sraeisid/FaNDOM/HCC827/HCC827_merged_110X_RawMolecules.bnx'
output_dir = '/nucleus/projects/sraeisid/FaNDOM/HCC827/raw_mol'
name = 'HCC827'
header = ''
print_line = ''
i = 0
counter= 1
with open(input_dir, 'r') as f:
    for line  in f :
        if line.startswith('#'):
            header = header + line
        else:
            i=i+1
            print_line = print_line + line
            if i == 1000000:
                with open(output_dir+'/'+name+'_'+str(counter)+'.bnx', 'w')as g:
                    g.write(print_line)
                    i=0
                    counter+=1
                    print_line = ''
    # print(header)
    with open(output_dir + '/' + name + '_' + str(counter) + '.bnx', 'w')as g:
        g.write(header)
        g.write(print_line)






