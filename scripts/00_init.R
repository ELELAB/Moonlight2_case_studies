# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Set working directory ---------------------------------------------------
setwd("./scripts")

# Check all req are install 
# Check packages ----------------------------------------------------------
# Checking for the required packages download if not present
install_if_needed(c("Moonlight2R",
                    "devtools",
                    "TCGAbiolinks"))

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

