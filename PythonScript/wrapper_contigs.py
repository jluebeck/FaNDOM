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
args = parser.parse_args()
c = 'hg19'
if args.chrom == '38':
    c = 'hg38'
else:
    c = 'hg19'

os.chdir(args.fandom)
out = args.output
if out[-1] != '/':
    out = out + '/'
name = args.name
pre_proces_cmd = 'python PythonScript/preprocess2.py -q ' + args.query + ' -o ' + out + args.name + '_preprocess'
os.system(pre_proces_cmd)
fandom_cmd = './FaNDOM -t=' + args.thread + ' -r=' + args.ref + ' -q=' + out + args.name + '_preprocess.cmap' + ' -sname=' + out + args.name + ' -outfmt=xmap '
os.system(fandom_cmd)
remove_part_cmd = 'python PythonScript/remove_part.py -p ' + out + args.name + '_partial.xmap -f ' + out + args.name + '.xmap -o ' + out + args.name + '_full'
os.system(remove_part_cmd)
filter_partial_cmd = 'python PythonScript/filter_contigs.py -i ' + out + args.name + '_partial.xmap -r ' + args.ref + ' -o ' + out + name + '_partial_filtered'
os.system(filter_partial_cmd)
filter_partial_post_cmd = 'python PythonScript/post_process.py -f ' + out + name + '_partial_filtered.xmap -d ' + out + args.name + '_preprocess_dic -o ' + out + name + '_partial_filtered'
os.system(filter_partial_post_cmd)
full_post_cmd = 'python PythonScript/post_process.py -f ' + out + args.name + '_full.xmap -d ' + out + args.name + '_preprocess_dic -o ' + out + args.name + '_full'
os.system(full_post_cmd)
assemble_cmd = 'python PythonScript/assemble_reads.py -i ' + out + args.name + '_full_post_process.xmap -o ' + out + args.name + '_full_post_process'
os.system(assemble_cmd)
os.chdir(out)
os.system(
    'cat ' + out + args.name + '_full_post_process_assembled.xmap ' + out + name + '_partial_filtered_post_process.xmap >' + name + '_all.xmap')
os.chdir( args.fandom)
assemble_cmd = 'python PythonScript/assemble_reads.py -i ' + out + name + '_all.xmap' + ' -o ' + out + name + '_all'
os.system(assemble_cmd)
sv_detect_cmd = 'python PythonScript/SV_detection_contigs.py -i ' + out + name + '_all_assembled.xmap -l 1 -c ' + c + ' -q ' + args.query + ' -o ' +out + 'SV'
os.system(sv_detect_cmd)
os.chdir( out)
os.system('cat '+ out + args.name + '_full_post_process.xmap '+ out + name + '_partial_filtered_post_process.xmap >' + name + '_all_indel.xmap')
os.chdir( args.fandom)
indel_detect_cmd = 'python3 PythonScript/indel_detection_contigs.py -r ' + args.ref + ' -c ' + c + ' -m ' + args.query + ' -a ' + out + name + '_all_indel.xmap -o ' + out + 'indel'
os.system(indel_detect_cmd)
os.chdir( out)
os.system("grep 'ins' indel > insertion")
os.system("grep 'del' indel > deletion")
# os.chdir( args.fandom)
# cluser_cmd = 'python PythonScript/cluster_indels.py -i ' + out +'insertion '+ '-o '+out+'insertion_clustered'
# os.system(cluser_cmd)
# cluser_cmd = 'python PythonScript/cluster_indels.py -i ' + out +'deletion '+ '-o '+out+'deletion_clustered'
# os.system(cluser_cmd)
