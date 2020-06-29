def coef(d , key):
    c = 3000
    key = int(key)
    id = d[key]
    ans = 0
    while(1):
        if key-1 in d.keys():
            if d[key-1] == id:
                key = key -1
                ans +=1
            else:
                break
        else:
            break
    return  c * ans

d = {}
with open('pre_process_dic','r') as f:
    for line in f :
        line = line.strip().split('\t')
        d[int(line[0])] = line[1]
with open('filtered.xmap','r') as f :
    with open('filtered_post_process.xmap', 'w') as g:
        for line in f :
            if line.startswith('#'):
                g.write(line)
            else:
                line2 = line.strip().split('\t')
                if int(line2[0]) in d.keys():
                    c = coef(d,line2[0])
                    line2[0] = d[int(line2[0])]
                    alignment = line2[-1]
                    alignment = alignment.split(')')[:-1]
                    new_align = ''
                    for i in alignment:
                        b = i.split(',')
                        new_align =new_align +  b[0]+','+str(int(b[1])+c)+')'
                    line2[-1] = new_align
                g.write('\t'.join(line2))
                g.write('\n')

