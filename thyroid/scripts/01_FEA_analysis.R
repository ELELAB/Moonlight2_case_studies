# This script makes the functional enrichment analysis and plot the results
# It should run prior to 'run_moonlight.R'

# SETUP ---------------------------------
# Library ---------------
library(Moonlight2R)
library(tidyverse)
# Data -----------------
DEA_Basal <- get(load("../data/rawdata/THCA_DEGs_NTvs.TP.rda"))
ind <- which(duplicated(DEA_Basal$ID))
DEA_Basal <- DEA_Basal[-ind,] %>% 
  rownames_to_column(var = "row_ID") %>% 
  column_to_rownames(var = "ID") 
DEA_Basal <- DEA_Basal[,c("logFC","AveExpr",
                          "t","P.Value","adj.P.Val","B")]

# Run functions -------------------------
FEA_Basal <- FEA(DEGsmatrix = DEA_Basal)
save(FEA_Basal, file = "../data/FEA_Thyroid.rda")


# Plot result ---------------------------
png(filename = "../results/plotFEA_Thyroid.png",
	width = 8, height = 4, units = "in", pointsize = 4, res = 1200)
plotFEA(FEA_Basal, topBP = 20)
dev.off()

