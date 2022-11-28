#This script will plot two heatmaps of the moonlight scores  

# SETUP ----------------------------------
# Library ---------------
library(Moonlight2R)

# Data ------------------
load("../data/URA_Basal.rda")
load("../results/Oncogenic_mediators_mutation_summary.rda")
load("../results/DEG_Mutations_Annotations.rda")

# PLOT 1 ---------------------------------
# Plot the most mutated driver genes 
plotMoonlight(DEG_Mutations_Annotations,
              Oncogenic_mediators_mutation_summary,
              dataURA = URA_Basal,
              gene_type = "drivers",
              additionalFilename = "../results/022_Drivers_")

# PLOT 2 ---------------------------------
# Plot the most mutated oncogenic mediators 
plotMoonlight(DEG_Mutations_Annotations,
              Oncogenic_mediators_mutation_summary,
              dataURA = URA_Basal,
              gene_type = "mediators",
              additionalFilename = "../results/022_Mediators_")
