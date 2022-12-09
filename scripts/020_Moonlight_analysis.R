# This script runs the moonlight pipeline

# SET UP -------------------------------
# Library ----------------
library(TCGAbiolinks)
library(Moonlight2R) 
library(tidyverse)

# Data -------------------
DEA_Basal <- get(load("../data/rawdata/BRCA_BRCA.Basal_dataDEGs.rda"))
Filt_Basal_raw <- get(load("../data/rawdata/BRCA_BRCA.Basal_dataFilt_HUGO.rda"))
Maf_Basal <- read.csv("../data/rawdata/mutations.csv")

# Wrangle data -------------------------
# Removing non-tumor data points: 
tum.samples <- TCGAquery_SampleTypes(barcode = colnames(Filt_Basal_raw), 
                                     typesample = "TP")
Filt_Basal <- Filt_Basal_raw[,which(colnames(Filt_Basal_raw) %in% tum.samples)]

print('everything is loaded')
# Call Moonlight functions -----------------------
#GRN -------
GRN_Basal <- GRN(TFs = rownames(DEA_Basal),
                 DEGsmatrix = DEA_Basal,
                 normCounts = Filt_Basal,
                 kNearest = 3,
                 DiffGenes = TRUE,
                 nGenesPerm = 1000,
                 nBoot = 100)
print('grn is finished')
save(GRN_Basal, file = "../data/GRN_Basal.rda")

#URA ------
URA_Basal <- URA(dataGRN = GRN_Basal,
                  DEGsmatrix = DEA_Basal,
                  nCores = 4,
                  BPname = c("proliferation of cells","apoptosis"))
save(URA_Basal, file = "../data/URA_Basal.rda")

# PRA ------
PRA_Basal <- PRA(dataURA = URA_Basal,
                  BPname = c("proliferation of cells", "apoptosis"),
                  thres.role = 0)
save(PRA_Basal, file = "../data/PRA_Basal.rda")

# DMA -------
Maf_Basal_filtered <- Maf_Basal %>%
  mutate(VAF = t_alt_count/t_depth) %>% 
  filter(VAF > 0.05, 
         t_alt_count >= 3,
         t_depth > 30,
         n_depth >10)

DMA_Basal <- DMA(dataMAF = Maf_Basal_filtered,
                 dataDEGs = DEA_Basal,
                 dataPRA = PRA_Basal,
                 results_folder = "../results",
                 coding_file = "/data/databases/CScape/CScape-20210624/css_coding.vcf.gz",
                 noncoding_file = "/data/databases/CScape/CScape-20210624/css_noncoding.vcf.gz")
save(DMA_Basal, file = "../data/DMA_Basal.rda")
