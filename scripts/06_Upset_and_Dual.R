# This script creates a Upset plot 
# The first plot is of oncogenic mediators 
# The second plot is of driver genes
# Based on these plots there is then made an enrichemnt analysis of all drivers
# and the possible dual-role driver genes


# Libraries -------------------------------
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(viridis)
library(Moonlight2R)
library(enrichR)
library(gridExtra)


source("99_functions.R")

# Load data -------------------------------
load("../results/Oncogenic_mediators_mutation_summary.rda")
#load("../results/DEG_Mutations_Annotations.rda")

# Mediator PLOT -----------------------------------------------------------------------------
## Wrangle data ----------------------------
names_sets <- c("Moonlight_OCG", "Moonlight_TSG", "NCG_OCG", "NCG_TSG", "NCG_Candidate")
Moonlight_OCG <- Oncogenic_mediators_mutation_summary %>% filter(Moonlight_Oncogenic_Mediator == "OCG") %>% pull(Hugo_Symbol)
Moonlight_TSG <- Oncogenic_mediators_mutation_summary %>% filter(Moonlight_Oncogenic_Mediator == "TSG") %>% pull(Hugo_Symbol)
NCG_OCG <- NCG %>% filter(NCG_driver == "OCG") %>% pull(symbol)
NCG_TSG <- NCG %>% filter(NCG_driver == "TSG") %>% pull(symbol)
NCG_Candidate <-  NCG %>% filter(NCG_driver == "Candidate") %>% pull(symbol)

myset1 <- list(Moonlight_OCG, Moonlight_TSG, NCG_OCG, NCG_TSG, NCG_Candidate)

## Upset Plot data -------------------------------
upset_plot(file_name = "06_upsetplot_oncogenic_mediators.pdf", dir_output = "../results/",
           sets_list = myset1, names_sets = names_sets, title_plot = "Oncogenic mediators",
           viridis_color = "viridis")


# Driver PLOT -----------------------------------------------------------------------------
## Wrangle data ----------------------------
Drivers <- Oncogenic_mediators_mutation_summary %>% filter( CScape_Driver >= 1)

names_sets <- c("Moonlight_OCG", "Moonlight_TSG", "NCG_OCG", "NCG_TSG", "NCG_Candidate")
Moonlight_OCG <- Drivers %>% filter(Moonlight_Oncogenic_Mediator == "OCG") %>% pull(Hugo_Symbol)
Moonlight_TSG <- Drivers %>% filter(Moonlight_Oncogenic_Mediator == "TSG") %>% pull(Hugo_Symbol)
NCG_OCG <- NCG %>% filter(NCG_driver == "OCG") %>% pull(symbol)
NCG_TSG <- NCG %>% filter(NCG_driver == "TSG") %>% pull(symbol)
NCG_Candidate <-  NCG %>% filter(NCG_driver == "Candidate") %>% pull(symbol)

myset2 <- list(Moonlight_OCG, Moonlight_TSG, NCG_OCG, NCG_TSG, NCG_Candidate)


## Plot data ----------------------------------------------------------------
upset_plot(file_name = "06_upsetplot_drivers.pdf", dir_output = "../results/",
           sets_list = myset2, names_sets = names_sets, title_plot = "Drivers Genes",
           viridis_color ="turbo")

## NCG vs Moonlight driver -------------------------------------------------
overlap_tab <- Drivers %>% filter(NCG_driver != 'Candidate') %>% 
  select(Moonlight_Oncogenic_Mediator, NCG_driver, Hugo_Symbol) %>% 
  group_by(Moonlight_Oncogenic_Mediator, NCG_driver) %>% 
  summarise(Genes = list(Hugo_Symbol)) 
write_csv(overlap_tab, file = "../results/06_overlap_moonlight_drivers_vs_ncg.csv")

## Driver Enrichment Analysis ----------------------------------------------
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", 
         "GO_Biological_Process_2021", "KEGG_2021_Human")

# General enrichment of driver genes
driver_genes <- Drivers %>% pull(Hugo_Symbol)
driver_enrich <- enrichr(driver_genes, databases = dbs)
p1 <- goplot(data = driver_enrich$GO_Molecular_Function_2021, title = "MF", top = 8)
p2 <- goplot(data = driver_enrich$GO_Cellular_Component_2021, title = "CC" , top =8)
p3 <- goplot(data = driver_enrich$GO_Biological_Process_2021, title = "BP", top =8)
p4 <- goplot(data = driver_enrich$ GO_Molecular_Function_2021, title = "KEGG", top =8)
grid.arrange(p1, p2, p3,p4, ncol=2, nrow = 2) #, top = textGrob("Molecular Function",gp=gpar(fontsize=20,font=3)))

# DUAL ROLE DRIVER GENES --------------------------------------------- 
dual_drivers <- append(Moonlight_TSG[Moonlight_TSG %in% NCG_OCG],
                       Moonlight_OCG[Moonlight_OCG %in% NCG_TSG])
x<- enrichr(genes = dual_drivers, databases = dbs)
y <- rbind(x[[1]], x[[2]],x[[3]], x[[4]]) # combine all three types

# Genes they agree on 
tum <- Moonlight_TSG[Moonlight_TSG %in% NCG_TSG]
onc <- Moonlight_OCG[Moonlight_OCG %in% NCG_OCG] 
onco_x <- enrichr(genes = onc, databases = dbs)
tum_x <- enrichr(genes=tum, databases = dbs)

#Molecular function
p1 <- goplot(data = x$GO_Molecular_Function_2021, title = "Dual role", top = 10)
p2 <- goplot(data = onco_x$GO_Molecular_Function_2021, title = "Oncogenes" , top =10)
p3 <- goplot(data = tum_x$GO_Molecular_Function_2021, title = "Tumor suppressors", top =10)
grid.arrange(p1, p2, p3,ncol=3, top = textGrob("Molecular Function",gp=gpar(fontsize=20,font=3)))

#Biological Process
p1 <- goplot(data = x$GO_Biological_Process_2021, title = "Dual role", top = 10)
p2 <- goplot(data = onco_x$GO_Biological_Process_2021, title = "Oncogenes" , top =10)
p3 <- goplot(data = tum_x$GO_Biological_Process_2021, title = "Tumor suppressors", top =10)
grid.arrange(p1, p2, p3,ncol=3, top = textGrob("Biological Process",gp=gpar(fontsize=20,font=3)))

#KEGG
p1 <- goplot(data = x$KEGG_2021_Human, title = "Dual role", top = 10)
p2 <- goplot(data = onco_x$KEGG_2021_Human, title = "Oncogenes" , top =10)
p3 <- goplot(data = tum_x$KEGG_2021_Human, title = "Tumor suppressors", top =10)
grid.arrange(p1, p2, p3,ncol=3) #, top = textGrob("KEGG Pathways",gp=gpar(fontsize=20,font=3)))




## Compare Dual drivers with Shen et al. 2018 -----------------
Shen <- read.delim(file = "../data/rawdata/Shen_dual_genes.txt",
                   header = FALSE)
inner_join(as.data.frame(dual_drivers), Shen, by = c('dual_drivers' = 'V1'))
inner_join(as.data.frame(Drivers), Shen, by =c('Hugo_Symbol' = 'V1')) 
