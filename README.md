
Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

# Moonlight2_DMA_basal_like

An Automatized Workflow to Study Mechanistic Indicators for Driver Gene Prediction with Moonlight
Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico, Catharina Olsen, Matteo Tiberti, Elena Papaleo*, submitted to biorxiv

contacts for the repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Astrid Saksager: abrisa-at-dtu.dk

PLEASE, CITE THE PUBLICATION ABOVE IF YOU USE THE SCRIPTS FOR YOUR OWN RESEARCH

This Project-folder contains the entire Basal-like Breast Cancer study
AIM: To find driver genes for Basal-like patients by running Moonlight on BRCA-Basal specific data from TCGA.

## Installing requirements and reproducing the analysis

All the analyses have been performed on a GNU/Linux server.

### CScape-somatic

In order to run the analyses you will need to have available the pre-calculated
[CScape-somatic](http://cscape-somatic.biocompute.org.uk) scores. These are
downloadable from CScape-somatic website. Please download at least the
`css_coding.vcf.gz` and `css_noncoding.vcf.gz` files and store them in a local
directory of your choice. Then modify the `Run_Basal.sh` script so that the
system variables defined therein refer to the location of these two files
(see the file itself for an example).

### Computing environment

In order to reproduce the paper data, you will need to set up a conda environment
on which the expected version of R and the required packages will be installed;
this requires being able to run Anaconda by means of the `conda` executable.

If you don't have access to `conda` please see the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html) on instructions on how to install Miniconda.

Once you have access to `conda`, you can

1. clone our github repository into a local directory on your local machine:

```
git clone https://github.com/ELELAB/Moonlight2_DMA_basal_like.git
cd Moonlight2_DMA_basal_like
```

2. create a virtual environment using conda activate it:

```
conda create --prefix ./env_Basal -c conda-forge r-base=4.2 r-pacman r-curl r-ragg r-renv r-osfr r-cairo
conda activate ./env_Basal
```

3. run the analysis:

```
$ bash ./Run_Basal.sh
```

**WARNING**: our script use the [renv](https://rstudio.github.io/renv/articles/renv.html)
R package to handle automatic dependency installation. `Renv` writes packages in
its own cache folder, which is by default in the user's home directory. This might
not be desireable if free space in the home directory is limited as this cache
takes about 2GB for this project. You can change the location of the Renv root
folder by setting a system environment variable - please see comments in the 
`Run_Basal.sh` script.

This script will:

1. Install in the environment all the necessary R packages, including Moonlight2R
from [our GitHub repository](https://www.github.com/ELELAB/Moonlight2R)
2. Download from [our OSF repository](https://osf.io/eq9wj/) the `data` folder
that contains data required for the analysis (see below). This takes approximately
700MB of disk space.
3. perform the analysis

## Content of the downloaded directories

### `./data`

This contains the data generated by Moonlight2:
`Basal_FEA.rda`, `Basal_URA.rda`, `Basal_GRN.rda`, `Basal_PRA.rda`, `Basal_DMA.rda`

### `./data/rawdata`

This data that is not generated in this project, but was used as input to
generate the results. It includes:
  - Mutation from TCGA-BRCA (`mutations.csv`)
  - Differential Gene expression from TCGA-BRCA (`BRCA_BRCA.Basal_dataDEGs.rda`)
  - Normal counts of gene expression (`BRCA_BRCA.Basal_dataFilt_HUGO.rda`)
  - Driver genes predicted by GUST (`GUST_BRCA.csv`)
  - List of transcription factors from TRRUST (`trrust_rawdata_human.tsv`)
  - List of human kinases from Kinhub (`kinhub.tsv`)
  - List of PAM50 genes (`pam50_list.txt`)
  - List of dual role driver genes from Shen et al. 2018 (`Shen_dual_genes.txt`)

### `./scripts/`
The script `00_init.R` downloads the raw data from our OSF respository, 
loads and installs the software packages required to run the rest of the scripts,
and performs all the necessary calculations by sourcing a number of other
scripts.

The script are numbered and the results from each are numbered accordingly.

### `./results/`
This folder contains plots and tables used in the paper. The numbers in filenames
refer to the script they are generated from.

## Paper contents

The paper contents corresponds to the following files:

| Figure/Table   | Filename                                              |
| -------------- | ----------------------------------------------------- |
| Figure 2       | `results/021_heatmap_complete.pdf`                    |
| Figure 3       | `results/022_Mediators_moonlight_heatmap.pdf`         |
| Figure 4       | `results/032_promoter_moonlight_heatmap.pdf`          |
| Figure 5       | `results/05_GO_KEGG_top_adjpval_OCG_TSG.png`          |
| Figure 6       | `results/06_upsetplot_drivers.pdf`                    |
| Table 1        | The CScape scores are derived from Rogers et al. 2020 |
| Table 2        | text output of `scripts/030_summarise.R`              |
| Table 3        | `results/032_promoter_mutation_genes.csv`             |
| Table 4        | `results/033_driver_TF.csv`                           |
| Table S2       | In OSF repository: `data/DMA_Basal.rda.gz`            |
| Table S3       | `results/06_overlap_moonlight_drivers_vs_ncg.csv`     |
| Table S4       | `results/07_overlapping_genes_moonlight_ncg_gust.csv` |
| Sup. Figure S1 | `results/05_variant_class_percent_driver_TSG_OCG.png` |
| Sup. Figure S2 | `results/98_mutations_per_patient.png`                |
