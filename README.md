# FaNDOM 
**Fa**st **N**ested **D**istance aligner for **O**ptical **M**aps

## About

FaNDOM performs alignment of Bionano Saphyr optical map molecules and contigs to a reference, using a seed-based filter.
FaNDOM is implemented in C++ and supports multithreading. Currently, FaNDOM is available for beta testing.

FaNDOM is developed by Siavash Raeisi Dehkordi and Jens Luebeck.

## Installation

FaNDOM requires `cmake 3.1` or higher. This should already be satisfied on most modern Unix systems. 
To check your cmake version, type `cmake --version`.  

To install FaNDOM, do

```
git clone https://github.com/jluebeck/FaNDOM
cd FaNDOM
cmake CMakeLists.txt && make 
```

## Usage

FaNDOM takes as input Bionano Saphyr molecules stored in `.bnx` or `.cmap` form or assembled contigs in `.cmap` form, and a pre-processed reference genome.
Right now, FaNDOM supports GRCh37(hg19) and GRCh38 available in `reference_genomes` folder.
FaNDOM outputs alignments of the OM molecules in FaNDOM's `.fda` or `.xmap` file format.

### Command line arguments

##### Required arguments
- `-r=` Path to reference genome for alignment
- `-q=` Path to Bionano Saphyr molecules or contigs
- `-sname=` Prefix for output files
- `-outfmt=` Specify output format for alignments. We support `.fda`and `.xmap` format.

##### Optional arguments (basic)
- `-multimap` Report multiple alignments per molecule (default: only highest scoring alignment)
- `-version` Print version and exit
- `-t=` Number of threads to use (recommend 12+)
- `-padding=` Additional size (in bp) around seed region to open alignment window (default: 1000)
- `-no_partial=` 0 if Just looking for full alignment (default: 1)
- `-ref38=` If you are using reference genome assembly GRCh38 set it to True (default: False)
- `-rescale=` True if data are raw molecule and it is necessary to rescale them (default: False)

##### Optional arguments (advanced)
- `-tolerance=` Seeding label position error tolerance (default: 350)
- `-rank=` Seed ranks to consider (default: 150)
- `-threshold=` Seed chain mininum length (default: 3)
- `-band_width=` Half of band width (default 6000)
### Example
`./FaNDOM -t=20 -r=reference_genomes/hg19_Merge_800_DLE.cmap -q=query/EXP_REFINEFINAL1_q.cmap -sname=output/output -outfmt=xmap`


For making sure that you install FaNDOM correctly, when you are at FaNDOM directory run the following command:
```
./FaNDOM -t=1 -r=test_data/reference.cmap -q=test_data/query.cmap -sname=test_data/res -outfmt=xmap

```


## Wrapper for SV analysis of assembled contig data
If you want to run whole the pipeline for detecting SVs on assembled contigs, you can use a python script in "Pythonscript" folder named `wrapper_contigs.py` 
-  `-q` Path to contigs file in '.cmap' file format.
-  `-r` Path to reference file. It should be in cmap format.
-  `-f` Path to foldar that contains executable file of FaNDOM.
-  `-t` Number of threads.
-  `-o` Path to a directory for saving all alignments and SV calls.
-  `-n` Name of alignment files
-  `-c` Assemble of reference that is used. 19 for GRCh37(hg19) and 38 for GRCh38(hg38)
-  `-m` If you are aligning contigs having more than 300 labels, use mode 1 to preprocess input data and generate shorter contigs, otherwise use mode 2. 
As an example:
```
python wrapper_contigs.py -f /home/FaNDOM -t 22 -r /home/reference_genomes/hg19_Merge_800_masked.cmap -q /home/exp_refineFinal1_sv/exp_refineFinal1_merged_q.cmap -n Covid_3_41 -o /home/res -c 19 -m 1
```
For making sure that you install FaNDOM correctly, when you are at FaNDOM directory run the following command:

```
python PythonScript/wrapper_contigs.py -f $PWD -t 1 -r test_data/reference.cmap -q test_data/query.cmap -n test -o $PWD/test_data -c 19 -m 1

```
It should run pipeline for simple datasets and have the results at `test_data/res` directory.
## Python scripts
The following scripts are used inside the SV wrapper - `wrapper_contigs.py`, and can be invoke separately if desired.
### `autorescale.py` script
This script used for finding the best rescale factor for molecules. We highly recommend to use this script if you are using raw molecules.
-  `-q` Path to molecules file. It can be bnx or cmap file.
-  `-r` Path to reference file. It should be in cmap format.
-  `-f` Path to foldar that contains executable file of FaNDOM.
-  `-t` Number of threads.
-  `-o` Path to a directory for saving the rescaled molecules file.
-  `-s` Down sample size. Default value is 500.
As an example:
```
python3.5 autorescale.py -q /data/molecule/a.bnx -r /data/ref/hg19_DLE.cmap -f /data/FaNDOM/ -t 10 -o /data/molecule/rescale
```
### `preprocess.py` script
For a case you want to align very large assembled contigs which can have more than 500 labels on one contigs, at first you need to run preprocess.py on your query to separate large contigs to smaller pices and then run Fandom on it.
-  `-q` Path to assembled contig file. It should be in '.cmap' format.
-  `-m` Maximum size for new separated molecules in terms of number of labels. Default value is 25.
-  `-o` Path to a directory for savind processed queries and also save a dictionary for converting Fandom alignment to actual molecules.

As an example:
```
python preprocess.py -q H460_DLE1_EXP_REFINEFINAL1.cmap -o /Output/processed2 -m 50
```
### `post_process.py` script
This script used for remap aligments to first molecule file. For doing that you need a file ending with 'dic' that preprocess.py script made.
-  `-f` Path to Fandom alignmrnt output.
-  `-d` Path to pre_process.py dictionary.
-  `-o` Output directory for translated alignments.
As an example:
```
python post_process.py -f Fandom_output.xmap -d processed_dic -o fandom_post_process
```
### `assemble_reads.py` script
This script is used for assembling partial alignments to a larger one( assembling each fragment of large contigs)
-  `-i` Path to unassembled alignmet file.
-  `-o` Output directory for assembled alignment file.
```
python assemble_reads.py -i Fandom_output.xmap -o Fandom_output_assembled.xmap
```

### `remove_part.py` script
This script is used for removing partial alignments from full alignments
-  `-f` Path to full alignments file in xmap format.
-  `-p` Path to partial alignments file in xmap format.
-  `-o` Output directory for only full alignments file.
```
python post_process.py -f Fandom_output.xmap -p Fandom_output_partial.xmap -o Fandom_output_full.xmap
```
### `filter_contigs.py` script
This script is used for filtering high confidence partial alignments for assembled contigs.
-  `-i` Path to partial alignments file in xmap format.
-  `-o` Output directory filtered partial alignments.
```
python filter_contigs.py -i Fandom_output_partial.xmap -o Fandom_output_partial_filtered.xmap
```
### `filter_individual.py` script
This script is used for filtering high confidence partial alignments for raw molecules.
-  `-r` Path to reference file. It should be in cmap format.
-  `-i` Path to partial alignments file in xmap format.
-  `-o` Output directory filtered partial alignments.
```
python filter_individual.py -r hg19_DLE.cmap -i Fandom_output_partial.xmap -o Fandom_output_partial_filtered.xmap
```
### `indel_detection_contigs.py` script
This script is used for finding indels in assembled contigs alignment files.
-  `-r` Path to reference file. It should be in cmap format.
-  `-a` Path to alignments file in xmap format.
-  `-o` Output directory for indel finding
-  `-m` Path to query file in bnx or cmap file format.
-  `-c` Assemble of reference that is used. 19 for GRCh37(hg19) and 38 for GRCh38(hg38)
-  `-g` Gene coordinate directory for human genome
```
python indel_detection_contigs.py -r hg19_DLE.cmap -g Gene_hg19.txt -a Fandom.xmap -o res/indel -c 19 -m query/query.contigs
```
### `indel_detection_individual.py` script
This script is used for finding indels in raw molecule alignment files.
-  `-r` Path to reference file. It should be in cmap format.
-  `-a` Path to alignments file in xmap format.
-  `-o` Output directory for indel finding
-  `-m` Path to query file in bnx or cmap file format.
-  `-c` Assemble of reference that is used. 19 for GRCh37(hg19) and 38 for GRCh38(hg38)
```
python indel_detection_individual.py -r hg19_DLE.cmap -a Fandom.xmap -o res/indel -c 19 -m query/query_contigs.cmap
```

### `SV_detection_contigs.py` script
This script used for detecting potential integration points.
-  `-i` Path to alignment file.
-  `-l` Minimum number molecules to support a integration point. Our suggestion for contigs is 1.
-  `-o` Output directory for list of integration points.
-  `-c` Assemble of reference that is used. 19 for GRCh37(hg19) and 38 for GRCh38(hg38)
-  `-q` Path to query file in cmap file format.
-  `-r` Path to reference file. It should be in cmap format.
-  `-g` Gene coordinate directory for human genome
As an example:
```
python SV_detection_contigs.py -i alignment.xmap -r hg19_DLE.cmap -g Gene_hg19.txt -l 1 -o SV.txt -q query/query_contigs.cmap -c 19
```

## `.fda` file format

We define a file format for FaNDOM output, where each alignment entry contains four lines (as defined in the header):

```
#0      ref_id  mol_id  aln_direction   ref_start_pos   ref_end_pos     mol_start_pos   mol_end_pos     mol_length
#1      total_score     mean_score      is_multimapped  is_secondary    aln_seed_num
#2      alignment [aln_index]:(ref_pos, mol_pos, mol_lab, score_delta)
#3      cigar [aln_index]:(delta_ref, delta_mol, mol_label_diff, delta_difference)
```

`is_secondary` is set to True if the molecule is multimapped (`is_multimapped = True`) and another alignmnet of the 
molecule in the file has a higher `total_score`. 

The `cigar` field specifies a list of tuples (tagged by the number in the alignment, starting at 0), with the following definitions:
- `delta_ref`: Distance in bp between this reference label in the alignment and the reference label in the alignment for the previous tuple.
- `delta_mol`: Distance in bp between this molecule label in the alignment and the molecule label in the alignment for the previous tuple. 
- `mol_label_diff`: Number of molecule labels between this cigar tuple and previous cigar tuple - 1 (number of skipped molecule labels).
- `delta_difference`: `delta_mol - delta_ref`.

