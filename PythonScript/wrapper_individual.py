import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fandom", help="FaNDOM dir", required=True)
parser.add_argument("-t", "--thread", help="Number thread", required=True)
parser.add_argument("-r", "--ref", help="Reference", required=True)
parser.add_argument("-q", "--query", help="Query", required=True)
parser.add_argument("-n", "--name", help="Name", required=True)
parser.add_argument("-o", "--output", help="Output_dir", required=True)
parser.add_argument("-c", "--chrom", help="Output_dir", required=True)
parser.add_argument("-m", "--minimum", help="Minimum Number of support", required=True)
args = parser.parse_args()
c = 'hg19'
gene_dir =''
if args.chrom == '38':
    c = 'hg38'
    gene_dir = args.fandom+'/reference_genomes/Gene_hg38.txt '
elif args.chrom == 'nh':
	c = 'nh'
else:
    c = 'hg19'
    gene_dir = args.fandom+'/reference_genomes/Gene_hg19.txt '
minimum_support = args.minimum

os.chdir(args.fandom)
out = args.output
if out[-1] != '/':
    out = out + '/'
name = args.name
os.chdir(out)
os.mkdir('molecules')
os.mkdir('alignments')
os.chdir( args.fandom)
print('Split molecules')
split_cmd = 'python PythonScript/split.py -o '+ out+'molecules/ -n '+ name + ' -i '+ args.query
os.system(split_cmd)
files = os.listdir(out+ 'molecules/')
# print(files)
print('Doing alignments')
for f in files:
	if f.startswith(name) and f.endswith('.bnx'):
		name_alignment = f[:-4]
		fandom_cmd = './FaNDOM -t=' + args.thread + ' -r=' + args.ref + ' -q=' + out + 'molecules/'+f + ' -sname=' + out + 'alignments/'+name_alignment +'_alignment'+ ' -outfmt=xmap -rescale'
		os.system(fandom_cmd)
os.chdir(out+'alignments')
os.system('cat *_partial.xmap > p.xmap')
os.system('cat *_alignment.xmap > f.xmap')
os.chdir(args.fandom)
print('Remove Partial alignments')
remove_part_cmd = 'python PythonScript/remove_part.py -p ' + out + '/alignments/p.xmap -f '+ out + 'alignments/f.xmap -o'  + out + '/alignments/'+name + '_full'
os.system(remove_part_cmd)
os.chdir(out+'alignments')
os.system(' cat p.xmap '+name+'_full.xmap >mix.xmap')
os.chdir(args.fandom)
print('Filter alignments')
filter_partial_cmd = 'python PythonScript/filter_individual.py -i ' + out + 'alignments/mix.xmap -o ' + out+'alignments/final_alignment_untranslate'
os.system(filter_partial_cmd)
print('SV detection')
sv_detect_cmd = 'python3 PythonScript/SV_detection_individual.py -i ' + out + 'alignments/final_alignment_untranslate.xmap -l '+minimum_support +' -c ' + c + ' -r=' + args.ref+ ' -q ' + args.query +' -g '+gene_dir+ ' -o ' +out + 'alignments/SV.txt'
if c=='nh':
	sv_detect_cmd = 'python3 PythonScript/SV_detection_individual.py -i ' + out + 'alignments/final_alignment_untranslate.xmap -l '+minimum_support +' -c ' + c + ' -r=' + args.ref+ ' -q ' + args.query + ' -o ' +out + 'alignments/SV.txt'
os.system(sv_detect_cmd)
translate_cmd = 'python3 PythonScript/translate.py -i '+ out+'alignments/final_alignment_untranslate.xmap -o '+ out +'alignments/final_alignment'
os.system(translate_cmd)
print('Done')

