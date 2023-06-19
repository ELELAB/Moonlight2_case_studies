#This script will plot two heatmaps of the moonlight scores  

# SETUP ----------------------------------
# Library ---------------
library(Moonlight2R)

# Data ------------------
URA_Thyroid <- get(load("../data/URA_Thyroid.rda"))
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))

# PLOT 1 ---------------------------------
# Plot the most mutated driver genes 
plotMoonlight(DEG_Mutations_Annotations,
              Oncogenic_mediators_mutation_summary,
              dataURA = URA_Thyroid,
              gene_type = "drivers",
              additionalFilename = "../results/022_Drivers_")

# PLOT 2 ---------------------------------
# Plot the most mutated oncogenic mediators 
plotMoonlight(DEG_Mutations_Annotations,
              Oncogenic_mediators_mutation_summary,
              dataURA = URA_Thyroid,
              gene_type = "mediators",
              additionalFilename = "../results/022_Mediators_")
