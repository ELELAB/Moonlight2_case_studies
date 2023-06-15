# This script runs the moonlight pipeline

# SET UP -------------------------------
# Library ----------------
library(TCGAbiolinks)
library(Moonlight2R) 
library(tidyverse)

options(bitmapType='cairo')


# Data -------------------
DEA_Basal <- get(load("../data/rawdata/LUAD_DEGs_NTvs.TP.rda"))
Filt_Basal_raw <- get(load("../data/rawdata/LUAD_dataFilt_HUGO.rda"))
Maf_Basal <- read.csv("../data/rawdata/mutations_LUAD.csv")

# Wrangle data -------------------------

# Wrangling DEA data to have genes in rows
ind <- which(duplicated(DEA_Basal$ID))
DEA_Basal <- DEA_Basal[-ind,] %>% 
  rownames_to_column(var = "row_ID") %>% 
  column_to_rownames(var = "ID") 
DEA_Basal <- DEA_Basal[,c("logFC","AveExpr",
                          "t","P.Value","adj.P.Val","B")]

# Removing non-tumor data points: 
tum.samples <- TCGAquery_SampleTypes(barcode = colnames(Filt_Basal_raw), 
                                     typesample = "TP")
Filt_Basal <- Filt_Basal_raw[,which(colnames(Filt_Basal_raw) %in% tum.samples)]

# Subsetting mutation data to have same samples as expression

# Obtain first three components of column names in filtered expression data
filt_basal_parts <- sub("^(\\w+-\\w+-\\w+).*", "\\1", colnames(Filt_Basal))

# Subset MAF based on matching patient barcodes with patient barcodes of expression barcodes
Maf_Basal <- Maf_Basal[grepl(paste(filt_basal_parts, collapse = "|"), Maf_Basal$Tumor_Sample_Barcode), ]

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
save(GRN_Basal, file = "../data/GRN_Lung.rda")

#URA ------
URA_Basal <- URA(dataGRN = GRN_Basal,
                  DEGsmatrix = DEA_Basal,
                  nCores = 4,
                  BPname = c("proliferation of cells","apoptosis"))
save(URA_Basal, file = "../data/URA_Lung.rda")

# PRA ------
PRA_Basal <- PRA(dataURA = URA_Basal,
                  BPname = c("proliferation of cells", "apoptosis"),
                  thres.role = 0)
save(PRA_Basal, file = "../data/PRA_Lung.rda")

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
                 coding_file = Sys.getenv('CSCAPE_CODING'),
                 noncoding_file = Sys.getenv('CSCAPE_NONCODING'))
save(DMA_Basal, file = "../data/DMA_Lung.rda")
