input_dir = '/nucleus/projects/sraeisid/Howard/de_novo_assembly/M249_DM/M249USR181005_SC1_mol_M249USR181005_SC1_RawMolecules.bnx'
output_dir = '/nucleus/projects/sraeisid/Howard/de_novo_assembly/M249_DM/raw_mol'
name = 'M249_DM'
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
    print(header)
    with open(output_dir + '/' + name + '_' + str(counter) + '.bnx', 'w')as g:
        g.write(header)
        g.write(print_line)






