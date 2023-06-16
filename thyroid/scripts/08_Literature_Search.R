# This script performs gene literature search of predicted driver genes

# Load libraries
library(Moonlight2R)
library(easyPubMed)
library(tidyverse)

# Load data
DMA_Thyroid <- get(load("../data/DMA_Thyroid.rda"))

# Wrangle data
drivers <- append(DMA_Thyroid$TSG, DMA_Thyroid$OCG)

# Perform literature search of driver genes
GLS_data <- GLS(genes = drivers)

# Save results
write_csv(GLS_data, file = "../results/08_driver_literature.csv") 
