# TRANSCRIPTION/PROMOTER MUTATIONS 

#Load libraries -------------------------
library(tidyverse)


# Load data ------------------------------
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
trrust <- read_tsv("./../data/rawdata/trrust_rawdata_human.tsv",col_names = FALSE) %>% 
  rename(TF = X1, Target = X2, Mode = X3, pumid =X4)


#Wrangle data -----------------------------
logfc <- DEG_Mutations_Annotations %>% select(Hugo_Symbol, logFC) %>% unique()
muts <- Oncogenic_mediators_mutation_summary %>% filter(CScape_Driver >= 1) %>% 
  select(Hugo_Symbol, CScape_Driver)
drivers <- left_join(muts, logfc)

trrust_tab <- trrust %>% filter(Mode == 'Activation'|
                                  Mode == 'Repression') %>% 
  inner_join(drivers, by = c("TF" = "Hugo_Symbol"))

trrust_tab %>% group_by(TF, Moonlight_Oncogenic_Mediator) %>% summarise(n =n())

# Save data ---------------------------------
write_csv(x = trrust_tab, 
          file = "../results/033_driver_TF.csv")
