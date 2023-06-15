# This script plots heatmaps of the DMA results 

# SETUP --------------------------------------
# Libraries --------------
library(Moonlight2R)
library(purrr)

# Data ------------------
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
PRA_Lung <- get(load("../data/PRA_Lung.rda"))

# Plotting -----------------------------------
# # PLOT 1 
# plotDMA(DEG_Mutations_Annotations,
#         Oncogenic_mediators_mutation_summary,
#         type = "complete",
#         additionalFilename = "../results/021_")
# 
# # PLOT 2 
# plotDMA(DEG_Mutations_Annotations,
#         Oncogenic_mediators_mutation_summary,
#         type = "split",
#         additionalFilename = "../results/")

# PLOT 3
Drivers <- Oncogenic_mediators_mutation_summary %>% 
  filter(CScape_Driver >= 1)

plotDMA(DEG_Mutations_Annotations,
        Oncogenic_mediators_mutation_summary = Drivers,
        type = "complete",
        additionalFilename = "../results/021_")

# PLOT 4 
# pam50 <-read.delim("../data/rawdata/pam50_list.txt")
# pam50 <- pam50$Basal 
# #load("../data/DMA_Basal.rda")
# #Oncogenic_Mediators <- append(DMA_Basal$TSG, DMA_Basal$OCG)
# #Oncogenic_Mediators <- names(append(PRA_Basal$TSG, PRA_Basal$OCG))
# plotDMA(DEG_Mutations_Annotations,
#         Oncogenic_mediators_mutation_summary,
#         genelist = pam50,
#         additionalFilename = "../results/Pam50_")
