# Additional Plots within the Moonlight package

# SETUP ---------------------------------------
# Libaries ---------------
library(Moonlight2R)

# Data -------------------
GRN_Lung <- get(load("../data/GRN_Lung.rda"))
URA_Lung <- get(load("../data/URA_Lung.rda"))
PRA_Lung <- get(load("../data/PRA_Lung.rda"))

# Plotting ------------------------------------
# GRN plot ------------
png(filename = "../results/plotGRN.png")
plotNetworkHive(GRN_Lung, knownDriverGenes, thres = 0.55)
dev.off()

# URA plot -------------
png(filename = "../results/plotURA.png")
plotURA(dataURA = URA_Lung[c(names(PRA_Lung$TSG), names(PRA_Lung$OCG)),, drop = FALSE],
        additionalFilename = "_apoptosis_proliferation")
dev.off()
