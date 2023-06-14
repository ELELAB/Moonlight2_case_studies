# This script performs gene literature search of predicted driver genes

# Load libraries
library(Moonlight2R)
library(easyPubMed)
library(tidyverse)

# Load data
DMA_Basal <- get(load("../data/DMA_Basal.rda"))

# Wrangle data
drivers <- append(DMA_Basal$TSG, DMA_Basal$OCG)

# Perform literature search of driver genes
GLS_data <- GLS(genes = drivers)

# Save results
write_csv(GLS_data, file = "../results/08_driver_literature.csv") 
