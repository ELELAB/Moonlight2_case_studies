
Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

# Moonlight2_case_studies

An Automatized Workflow to Study Mechanistic Indicators for Driver Gene Prediction with Moonlight
Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico, Catharina Olsen, Matteo Tiberti, Elena Papaleo*, 
bioRxiv 2022.11.18.517066; doi: https://doi.org/10.1101/2022.11.18.517066

contacts for the repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Astrid Saksager: abrisa-at-dtu.dk

PLEASE, CITE THE PUBLICATION ABOVE IF YOU USE THE SCRIPTS FOR YOUR OWN RESEARCH

This repository contains case studies associated with the new release of the MoonlightR package, namely Moonlight2R,
which implements new features associated with mutational analysis in a cancer cohort to find driver genes. The case
studies are conducted on Basal-like breast cancer, lung cancer, and thyroid cancer using data from The Cancer Genome
Atlas (TCGA). The scripts and results associated with these case studies are stored in this repository.

AIM: To find driver genes for Basal-like, lung, and thyroid cancer patients by running Moonlight2 on data from TCGA.

Below are instructions on the first steps to reproducing the analyses. Afterwards, please see README files in each 
cancer folder for specifics on how to reproduce analyses for the given cancer type.

## Installing requirements and reproducing the analysis

All the analyses have been performed on a GNU/Linux server.

NB: When reproducing the analyses and results, the user cannot expect to obtain identical results to the ones
in this GitHub repository and associated with the publication due to stochastic processes in the GRN step of the Moonlight protocol. 

### Computing environment

In order to reproduce the paper data, you will need to set up a conda environment
on which the expected version of R and the required packages will be installed;
this requires being able to run Anaconda by means of the `conda` executable.

If you don't have access to `conda` please see the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html) for instructions on how to install Miniconda.

Once you have access to `conda`, you can

1. clone our github repository into a local directory on your local machine:

```
git clone https://github.com/ELELAB/Moonlight2_case_studies.git
cd Moonlight2_case_studies
```

2. create a virtual environment using conda and activate it. 
The environment directory should be placed in the Moonlight2_case_studies folder:

```
conda env create --prefix ./env_Moonlight --file conda_env.yaml 
conda activate ./env_Moonlight
```

3. follow instructions in README of specific cancer type to reproduce analyses
associated with that cancer type

