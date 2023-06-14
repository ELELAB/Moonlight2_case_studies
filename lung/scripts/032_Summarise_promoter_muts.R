# TRANSCRIPTION/PROMOTER MUTATIONS 

#Load libraries -------------------------
library(Moonlight2R)
library(tidyverse)
library(ComplexHeatmap)
library(tidyHeatmap)
library(circlize)

# Load data ------------------------------
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
URA_Lung <- get(load("../data/URA_Lung.rda"))

# Load functions --------------------------
source("99_functions.R")

#Wrangle data -----------------------------
promoter_genes_tab <- DEG_Mutations_Annotations %>% filter(!is.na(Moonlight_Oncogenic_Mediator), 
                                                           CScape_Mut_Class == "Driver",
                                                           Potential_Effect_on_Transcription == 1,
                                                           Annotation == 'Promoter') %>% 
  select(ID, Chromosome, Start_Position, End_Position, Hugo_Symbol, Moonlight_Oncogenic_Mediator, logFC, 
         Variant_Type, Variant_Classification, CScape_Coding_score, CScape_Noncoding_score, 
         t_depth, n_depth, Annotation_Start, Annotation_End )

promoter_genes <-promoter_genes_tab %>% pull(Hugo_Symbol)


# Plot data ---------------------------------
plotMoonlight2(DEG_Mutations_Annotations,Oncogenic_mediators_mutation_summary,
               URA_Lung, genelist = promoter_genes,
               additionalFilename = '../results/032_promoter_')
# Save data ---------------------------------
write_csv(x = promoter_genes_tab, 
          file = "../results/032_promoter_mutation_genes.csv")
