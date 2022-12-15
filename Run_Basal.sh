#!/bin/bash

# uncomment and modify as necessary to set the location of the renv
# root directory. Default is $HOME/.cache/R/renv
# export RENV_PATHS_ROOT=/data/user/teo/scratch/renv

# Location of the CScape files - to be changed according to your system
export CSCAPE_CODING="/data/databases/CScape/CScape-20210624/css_coding.vcf.gz"
export CSCAPE_NONCODING="/data/databases/CScape/CScape-20210624/css_noncoding.vcf.gz"

Rscript scripts/00_init.R
