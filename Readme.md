## NMD-leafcutter

### Repository format
```
.
├── config/
│   ├── AT2.txt
│   ├── DD.txt
│   ├── FB.txt
│   └── MV.txt
├── data
│   └── bams
├── downstream_analysis/
│   ├── leafcutter_pipeline_report_donor.Rmd
│   └── leafcutter_top100_genes_report.Rmd
├── leafcutter/
├── logs/
├── REF/
│   ├── get_exon.sh
│   ├── get_reference.sh
│   ├── Readme.md
├── renv/
├── results/
├── src/
│   ├── run_leafcutter.py
│   ├── run_leafcutter_donor_covariate.py
│   ├── run_leafcutter_donor_covariate.sbatch
│   ├── run_all_celltypes.sh
│   └── compare_leafcutter_runs.py
├── Makefile
├── Readme.md
├── renv.lock
└── splicing_leafcutter.yaml

```
### About
This repository is designed to facilitate the analysis of Nonsense-Mediated Decay (NMD) BAM files (of cell types AT2, DD, Fibroblasts (FB) and Microvascular cells (MV)), focusing on differential splicing across cell type.

### Features
* Leafcutter Analysis: Performs differential splicing using leafcutter/scripts/leafcutter_ds.R
* Donor covariate: `run_leafcutter_donor_covariate.py` adds donor as a covariate in leafcutter_ds.R for the paired DMSO/Smg1i design
* Input Support: Handles direct BAM inputs
* Ouput Support: Stores the results in output_path of each cell run
* SLURM support: runs all cell types in parallel as a job array
* Reports: R Markdown reports for the pipeline audit and the top 100 significant genes

### Requirements:
* R v4.2.3
* conda
* samtools
* regtools
* rstan 2.21 and StanHeaders 2.21 (newer versions do not compile the leafcutter R package)

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

__5. install the leafcutter R package__
`leafcutter_ds.R` needs the leafcutter R package, which is built from the `leafcutter/leafcutter` folder of the clone.

Install the dependencies:
```{sh}
conda install -c conda-forge "r-rstan=2.21" "r-stanheaders=2.21" r-rcppparallel r-rstantools \
  r-shiny r-shinyjs r-dt r-shinycssloaders r-optparse r-intervals r-domc r-reshape2 r-gridextra r-r.utils
conda install -c conda-forge -c bioconda bioconductor-biobase bioconductor-dirichletmultinomial
Rscript -e 'install.packages(c("oompaBase", "oompaData", "TailRank"), repos="https://packagemanager.posit.co/cran/2024-03-01")'
```

Add RcppParallel to the package build (needed for the TBB headers used by Stan):
```{sh}
cd leafcutter/leafcutter
sed -i 's/^LinkingTo:/LinkingTo: RcppParallel (>= 5.0.1),/' DESCRIPTION
cat >> src/Makevars <<'EOF'

PKG_CPPFLAGS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "RcppParallel::CxxFlags()")
PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "RcppParallel::RcppParallelLibs()") $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "StanHeaders:::LdFlags()")
EOF
R CMD INSTALL .
cd ../..
```

Check the install:
```{sh}
Rscript -e 'library(leafcutter); cat("ok\n")'
```
If loading fails with a missing `libicu*.so` error, reinstall stringi with `conda install -c conda-forge r-stringi --force-reinstall`.

### Run the leafcutter script
Output can be pointed to the ```results``` folder. Sample run command:

```{sh}
python src/run_leafcutter.py --bam_file_list config/cell_type.txt --output_dir results/cell_type --cell_type --exon_file REF/exon_file.txt.gz
```

### Run the leafcutter script with donor as a covariate
Run from the repository folder. All paths can be relative, and the Leafcutter scripts are found automatically from the `leafcutter/` folder.

One cell type:
```{sh}
python src/run_leafcutter_donor_covariate.py \
  --bam_file_list config/FB.txt \
  --output_dir leafcutter_donor_results_v1/FB \
  --cell_type FB \
  --exon_file REF/gencode.v45.primary_assembly.txt.gz \
  --threads 8
```

All cell types, one after another (inside an `srun` session, with the conda environment active):
```{sh}
THREADS=8 bash src/run_all_celltypes.sh
```

All cell types in parallel on SLURM (one array task per file in `config/`):
```{sh}
mkdir -p logs
N=$(ls config/*.txt | wc -l)
sbatch --array=0-$((N-1)) src/run_leafcutter_donor_covariate.sbatch
```
Submit from the repository folder. Logs are written to `logs/`.

Both runner scripts process every `.txt` file in `config/`, and take the cell type name from the file name (`FB.txt` → `FB`).

### Data file path
* bam files are expected to be kept in ```data/bams```

### Information regarding input arguments
* ```--bam_file_list:``` should be the path to bam text files containing the bam paths and respective control vs treatment
* ```--output_dir:``` should be path to the output directory
* ```--cell_type :``` should be the cell type to be processed (AT2, DD, Fibroblasts (FB) and Microvascular (MV) cells)
* ```--exon_file:``` should be the path to the exon file generated 

Additional arguments for `run_leafcutter_donor_covariate.py`:
* ```--strand:``` regtools strand setting, `RF` (default) for fr-firststrand libraries
* ```--threads:``` threads for leafcutter_ds.R (defaults to `SLURM_CPUS_PER_TASK`, or 1)
* ```--min_samples_per_intron:``` minimum samples per intron for leafcutter_ds.R (default 3)
* ```--leafcutter_dir:``` path to the leafcutter clone, if it is not in the repository folder
* ```--cluster_script``` / ```--ds_script:``` override a single Leafcutter script, for example a patched DS script

### Preparing the BAM File List for Differential Splicing (Smg1i vs DMSO)

* When setting up a comparison between two conditions (e.g., Smg1i and DMSO), create a tab-separated text file (bam_groups.tsv) where:

* Each line contains the absolute path to a .bam file

* The second column specifies the experimental condition

* The first group (Smg1i) must appear before the second group (DMSO)

* For `run_leafcutter_donor_covariate.py`, a third column gives the donor, which is used as a covariate

Example file: ```bam_groups.tsv```

```
sample1.bam	Smg1i
sample2.bam	Smg1i
sample1.bam	DMSO
sample2.bam	DMSO

```

Example with donor: ```config/FB.txt```

```
data/bams/FB001V_Smg1i.markdup.sorted.bam	Smg1i	001V
data/bams/FB027U_Smg1i.markdup.sorted.bam	Smg1i	027U
data/bams/FB001V_DMSO.markdup.sorted.bam	DMSO	001V
data/bams/FB027U_DMSO.markdup.sorted.bam	DMSO	027U
```

### Important Notes

Order matters: Always list the first condition (Smg1i) before the reference/control condition (DMSO).

The file should be tab-delimited, not space-delimited.

File paths must be absolute for `run_leafcutter.py`. `run_leafcutter_donor_covariate.py` also accepts paths relative to the folder it is run from.

Ensure BAM files are sorted and indexed (each should have a corresponding .bam.bai file).

Sign of deltapsi: leafcutter_ds.R treats the first group in the file as the reference, so with Smg1i listed first, `deltapsi` is DMSO − Smg1i. A negative value means the intron is used more in Smg1i. The reports convert this to Smg1i − DMSO.

### Output files
the differential splicing analysis results are stored in two files
1. ```results/cell_type/cell_type_cluster_significance.txt```
2. ```results/cell_type/cell_type_effect_sizes.txt```

The donor covariate version also writes:
* ```junction_files/*.junc```: regtools junctions per sample (reused on a rerun; delete the folder if regtools settings change)
* ```samples.tsv```: groups file passed to leafcutter_ds.R (sample, treatment, donor)
* ```cell_type_perind_numers.counts.gz```: intron counts per sample
* ```logs/cell_type.log```: pipeline log

### Reports
Run from ```downstream_analysis/```:

```{sh}
Rscript -e 'rmarkdown::render("leafcutter_pipeline_report_donor.Rmd", params = list(results_dir = "../leafcutter_donor_results_v1"))'
Rscript -e 'rmarkdown::render("leafcutter_top100_genes_report.Rmd", params = list(results_dir = "../leafcutter_donor_results_v1"))'
```
The top 100 report also writes CSV tables to ```top100_tables/```.

### Comparing runs with and without the donor covariate
```{sh}
python src/compare_leafcutter_runs.py \
  --old_dir results \
  --new_dir leafcutter_donor_results_v1 \
  --cell_types AT2 DD FB MV \
  --min_dpsi 0.1 \
  --out comparison_out
```