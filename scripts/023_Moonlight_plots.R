# Additional Plots within the Moonlight package

# SETUP ---------------------------------------
# Libaries ---------------
library(MoonlightR)
# Data -------------------
load("../data/GRN_Basal.rda")
load("../data/URA_Basal.rda")
load("../data/PRA_Basal.rda")

# Plotting ------------------------------------
# GRN plot ------------
png(filename = "../results/plotGRN.png")
plotNetworkHive(GRN_Basal, knownDriverGenes, thres = 0.55)
dev.off()

# URA plot -------------
png(filename = "../results/plotURA.png")
plotURA(dataURA = URA_Basal[c(names(PRA_Balsa$TSG), names(PRA_Basal$OCG)),, drop = FALSE], 
        additionalFilename = "_apoptosis_proliferation")
dev.off()