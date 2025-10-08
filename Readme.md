## NMD-leafcutter

### Repository format
```
.
├── config/
│   ├── FB.txt
│   └── MV.txt
├── data
│   └── bams
├── REF/
│   ├── get_exon.sh
│   ├── get_reference.sh
│   ├── Readme.md
├── renv/
├── results/
├── src/
│   └── run_leafcutter.py
├── Makefile
├── Readme.md
├── renv.lock
└── splicing_leafcutter.yaml

```
### About
This repository is designed to facilitate the analysis of Nonsense-Mediated Decay (NMD) BAM files (of cell types, Fibroblasts and Microvascular cells), focusing on differential splicing across cell type.

### Features
* Leafcutter Analysis: Performs differential splicing using leafcutter/scripts/leafcutter_ds.R
* Input Support: Handles direct BAM inputs
* Ouput Support: Stores the results in output_path of each cell run

### Requirements:
* R v4.2.3
* conda

### Installation
__1. Clone the repository__ 

```{sh}
git clone git@github.com:NNeuralDynamics/nmd-leafcutter.git
cd nmd-leafcutter
```
__2. Run make to create a conda and r environment__ 
```{sh}
make setup
```
```{sh}
conda activate splicing
```
```
make install_renv
```
__3. get the GTF__
To get the gtf navigate to the ```REF``` folder and run ```get_reference.sh``` with the argument corresponding to the version number

```{sh}
cd REF
sh get_reference.sh 45
```
where 45 is the GTF version

__4. generate the exon_file.txt.gz__
To get the exon navigate to the ```REF``` folder and run ```get_exon.sh``` with the argument corresponding to the gtf file and output exon file

```{sh}
cd REF
sh get_exon.sh input.gtf.gz output.txt.gz
```
### Run the leafcutter script
Output can be pointed to the ```results``` folder. Sample run command:

```{sh}
python src/run_leafcutter.py --bam_file_list config/cell_type.txt --output_dir results/cell_type --cell_type --exon_file REF/exon_file.txt.gz
```
### Data file path
* bam files are expected to be kept in ```data/bams```

### Information regarding input arguments
* ```--bam_file_list:``` should be the path to bam text files containing the bam paths and respective control vs treatment
* ```--output_dir:``` should be path to the output directory
* ```--cell_type :``` should be the cell type to be processed (Fibroblasts (FB) and Microvascular (MV) cells)
* ```--exon_file:``` should be the path to the exon file generated 

### Preparing the BAM File List for Differential Splicing (Smg1i vs DMSO)

* When setting up a comparison between two conditions (e.g., Smg1i and DMSO), create a tab-separated text file (bam_groups.tsv) where:

* Each line contains the absolute path to a .bam file

* The second column specifies the experimental condition

* The first group (Smg1i) must appear before the second group (DMSO)

Example file: ```bam_groups.tsv```

```
sample1.bam	Smg1i
sample2.bam	Smg1i
sample1.bam	DMSO
sample2.bam	DMSO

```
### Important Notes

Order matters: Always list the first condition (Smg1i) before the reference/control condition (DMSO).

The file should be tab-delimited, not space-delimited.

File paths must be absolute, not relative.

Ensure BAM files are sorted and indexed (each should have a corresponding .bam.bai file).