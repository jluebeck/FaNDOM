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
-  `-w=` Word size of seed (default: 3)
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