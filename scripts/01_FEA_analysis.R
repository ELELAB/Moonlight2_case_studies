# This script makes the functional enrichment analysis and plot the results
# It should run prior to 'run_moonlight.R'

# SETUP ---------------------------------
# Library ---------------
library(Moonlight2R)
# Data -----------------
DEA_Basal <- get(load("../data/rawdata/BRCA_BRCA.Basal_dataDEGs.rda"))

# Run functions -------------------------
FEA_Basal <- FEA(DEGsmatrix = DEA_Basal)
save(FEA_Basal, file = "../data/FEA_Basal.rda")


# Plot result ---------------------------
png(filename = "../results/plotFEA_Basal.png",
	width = 8, height = 4, units = "in", pointsize = 4, res = 1200)
plotFEA(FEA_Basal, topBP = 20)
dev.off()

