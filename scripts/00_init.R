# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Set working directory ---------------------------------------------------
#setwd()

# Load repository 
devtools::install_github(repo = "ELELAB/MoonlightR2@feature_docs_update", 
                         auth_token = "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")

# Check all req are install 
# Check packages ----------------------------------------------------------
# Checking for the required packages download if not present
install_if_needed(c("MoonlightR",
                    "devtools",
                    "TCGAbiolinks"))

# Run Scripts
source("run_FEA.R")
source("plot_FEA.R")

source("run_moonlight.R")
#source("metadata_analysis.R") 
#source("minor_calculations.R")

source("plot_moonlight.R")
source("plot_DMA.R")
#source("plots.R")
#source("plot_additional.R")
source("plot_upset.R") 

