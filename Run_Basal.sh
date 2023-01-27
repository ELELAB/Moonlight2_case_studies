#!/bin/bash

# uncomment and modify as necessary to set the location of the renv
# root directory. Default is $HOME/.cache/R/renv
# export RENV_PATHS_ROOT=/data/user/teo/scratch/renv

# Location of the CScape files - to be changed if necessary.
# if a relative path is specified, it must be relative to the
# scripts folder inside this repository

export CSCAPE_CODING="$(pwd)/css_coding.vcf.gz"
export CSCAPE_NONCODING="$(pwd)/css_noncoding.vcf.gz"

Rscript scripts/00_init.R
