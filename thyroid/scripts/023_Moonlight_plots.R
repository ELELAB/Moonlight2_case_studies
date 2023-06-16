# Additional Plots within the Moonlight package

# SETUP ---------------------------------------
# Libaries ---------------
library(Moonlight2R)

# Data -------------------
GRN_Thyroid <- get(load("../data/GRN_Thyroid.rda"))
URA_Thyroid <- get(load("../data/URA_Thyroid.rda"))
PRA_Thyroid <- get(load("../data/PRA_Thyroid.rda"))

# Plotting ------------------------------------
# GRN plot ------------
png(filename = "../results/plotGRN.png")
plotNetworkHive(GRN_Thyroid, knownDriverGenes, thres = 0.55)
dev.off()

# URA plot -------------
png(filename = "../results/plotURA.png")
plotURA(dataURA = URA_Thyroid[c(names(PRA_Thyroid$TSG), names(PRA_Thyroid$OCG)),, drop = FALSE],
        additionalFilename = "_apoptosis_proliferation")
dev.off()
