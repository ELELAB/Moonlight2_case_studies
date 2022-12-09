library(pacman)

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Check all req are install 
# Check packages ----------------------------------------------------------
# Checking for the required packages download if not present

p_load('BiocManager')
p_load('osfr')
p_load('maftools')
p_load('enrichR')
p_load('UpSetR')
p_load('liftOver')
p_load_gh('ELELAB/Moonlight2R')
p_load('ggplot2')
p_load('tidyverse')

# Download data from OSF --------------------------------------------------
node = osf_retrieve_node('eq9wj')
files = osf_ls_files(node)
osf_download(files, recurse=TRUE, conflicts="skip")

# Set working directory ---------------------------------------------------
setwd("./scripts")

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

## Supplementary text 2:
source("98_MAF_quality_analysis.R")

