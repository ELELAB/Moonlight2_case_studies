# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Check all req are installed 
# Check packages ----------------------------------------------------------

# Download data from OSF --------------------------------------------------
library('osfr')
node = osf_retrieve_node('eq9wj')
files = osf_ls_files(node)
osf_download(files, recurse=TRUE, conflicts="skip")

# Set working directory ---------------------------------------------------
setwd("./scripts")

# restore packages from lockfile
library(renv)
renv::settings$external.libraries(c('../../env_Moonlight/lib/R/library/'))
renv::restore()
renv::activate()
detach('package:renv', unload=TRUE)

# load required packages
library('maftools')
library('enrichR')
library('UpSetR')
library('liftOver')
library('Moonlight2R')
library('ggplot2')
library('tidyverse')
library('easyPubMed')

# set up for Cairo rendering, for headless machines
options(bitmapType='cairo')

#Run scripts

## Moonlight
source("01_FEA_analysis.R")
source("020_Moonlight_analysis.R")
source("021_Moonlight_plotDMA.R")
source("022_Moonlight_plotMoonlight.R")
source("023_Moonlight_plots.R")

## Exploring specific genes and mutations
source("030_Summarise.R")
source("031_Summarise_plots.R")
source("032_Summarise_promoter_muts.R")
source("033_Summarise_trrust.R")

## Clincal
source("04_Metadata_analysis.R")

## Similarities and differences between TSG and OCG
source("05_OCG_vs_TSG.R")

## Comparing driver genes to other findings
source("06_Upset_and_Dual.R")
source("07_Breast_Cancer.R")

## Literature search of driver genes
source("08_Literature_Search.R")

## Supplementary text 2:
source("98_MAF_quality_analysis.R")

