# FaNDOM 
**Fa**st **N**ested **D**istance aligner for **O**ptical **M**aps

### About

FaNDOM performs alignment of Bionano Saphyr optical map molecules to a reference, using a seed-based filter.
FaNDOM is implemented in C++ and supports multithreading. Currently, FaNDOM is available for beta testing.

FaNDOM is developed by Siavash Raisi and Jens Luebeck.

### Installation

FaNDOM requires `cmake 3.1` or higher. This should already be satisfied on most modern Unix systems. 
To check your cmake version, type `cmake --version`.  

To install FaNDOM, do

```
git clone https://github.com/jluebeck/FaNDOM
cd FaNDOM
cmake && make 
```

### Usage

FaNDOM takes as input Bionano Saphyr molecules stored in `.bnx` or `.cmap` form, and a pre-processed reference genome.

FaNDOM outputs alignments of the OM molecules in FaNDOM's `.fda` file format. `.xmap` support coming soon.

#### Command line arguments

##### Required arguments
- `-r=` Path to reference genome for alignment
- `-q=` Path to Bionano Saphyr molecules
- `-s=` Prefix for output files

##### Optional arguments (basic)
- `-multimap` Report multiple alignments per molecule (default: only highest scoring alignment)
- `-version` Print version and exit
- `-t=` Number of threads to use (recommend 12+)
- `-padding=` Additional size (in bp) around seed region to open alignment window (default: 1000)

##### Optional arguments (advanced)
- `-tolerance=` Seeding label position error tolerance (default: 350)
- `-rank=` Seed ranks to consider (default: 300)
- `-threshold=` Seed chain mininum length (default: 3)

### `.fda` file format

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
-  `-q` Path to assembled contig file. It should be in cmap format.
-  `-m` Maximum size for new separated molecules in terms of number of labels. Default value is 100.
-  `-o` Path to a directory for savind processed queries and also save a dictionary for converting Fandom alignment to actual molecules.

As an example:
```
python preprocess.py -q H460_DLE1_EXP_REFINEFINAL1.cmap -o /Output/processed2 -m 200
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
### `SV_detect.py` script
This script used for detecting potential integration points.
-  `-a` Path to alignment output.
-  `-m` Minimum number molecules to support a integration point.
-  `-o` Output directory for list of integration points.
As an example:
```
python SV_detect.py -a alignment.xmap -m 1 -o SV.txt
```
