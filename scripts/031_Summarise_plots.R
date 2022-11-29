# This script generates plots of different parts of the results from DMA data 
# including mutations per patients, genes of interest ....


# Load libraries ------------------------------
library(tidyverse)
library(maftools)

# Load data -----------------------------------
load("../results/Oncogenic_mediators_mutation_summary.rda")
load("../results/DEG_Mutations_Annotations.rda")


# PLOTS -----------------------------------------------------------------------
# Mutations per Gene with CScape facet 
Oncogenic_mediators_mutation_summary %>% 
    pivot_longer(cols = c(CScape_Driver, CScape_Passenger, CScape_Unclassified),
                 names_to = "Cscape", values_to = "numbers") %>% 
  ggplot(aes(y = Cscape, x = numbers, fill = Moonlight_Oncogenic_Mediator), alpha = 0.8) +
  geom_boxplot()+
  theme_bw()

cscape_muts_om<-Oncogenic_mediators_mutation_summary %>% 
  pivot_longer(cols = c(CScape_Driver, CScape_Passenger, CScape_Unclassified),
               names_to = "Cscape", values_to = "Mutations_Per_Gene") %>% 
  ggplot(aes(x = Mutations_Per_Gene, fill = Moonlight_Oncogenic_Mediator), alpha = 0.8) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  facet_grid(cols = vars(Cscape))+
  theme_bw()
ggsave(cscape_muts_om, filename = "../results/031_cscape_mutations_per_OM.png")

# Numbers of mutations per patient -------------------

## Boxplot with colored outliers: 
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.95) + 1.5 * IQR(x)) #0.75 to get 8 patients
}

x <- DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
                                       into = "Patient_ID", sep = 12) %>%
  group_by(Patient_ID, CScape_Mut_Class) %>% 
  summarise(Mutations_per_Patient=n()) %>% 
  filter(CScape_Mut_Class != 'No_mutations') %>% 
  group_by(CScape_Mut_Class) %>% 
  mutate(outlier = ifelse(is_outlier(Mutations_per_Patient), 
                          Patient_ID, 
                          as.numeric(NA))) 
outs <-x %>%  filter(!is.na(outlier))

muts_per_patient_cscape <- x %>% 
  ggplot(aes(x = Mutations_per_Patient, y = CScape_Mut_Class)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_point(data = outs, aes(color = Patient_ID, size =3), show.legend = c(color =TRUE,size = FALSE))+
  #geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3, vjust = "top", angle = 45)+
 # xlim(0,420)+
  scale_x_continuous(limits = c(0,420), breaks = round(seq(0, 420, by = 20),1)) + 
  theme_bw()
ggsave(muts_per_patient_cscape, filename =  "../results/031_boxplot_mutations_per_patient_cscape.png")

## Boxplot total mutations in patients ----
y<-DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
                                       into = "Patient_ID", sep = 12) %>%
  mutate(Total = "Total") %>% filter(CScape_Mut_Class != 'No_mutations') %>%
  group_by(Patient_ID, Total) %>% summarise(Mutations_per_Patient=n()) %>% 
  group_by(Total) %>% 
  mutate(outlier = ifelse(is_outlier(Mutations_per_Patient), 
                          Patient_ID, 
                          as.numeric(NA))) 
outs_y <- y %>%  filter(!is.na(outlier))

muts_per_patient_total <- y %>% 
  ggplot(aes(x = Mutations_per_Patient, y = Total)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_point(data = outs_y, aes(color = Patient_ID, size =3), 
             show.legend = c(color =TRUE,size = FALSE))+
  scale_x_continuous( breaks = round(seq(0, 1000, by = 40),1)) + 
  theme_bw()
ggsave(muts_per_patient_total, filename =  "../results/031_boxplot_mutations_per_patient_total.png")


#  Bar plot per patient of their  number mutations
muts_per_patient_bar <- DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
                                         into = "Patient_ID", sep = 12) %>%
    drop_na(Patient_ID) %>% 
    ggplot(aes(x= Patient_ID, fill = CScape_Mut_Class), alpha = 0.8) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme_bw()
ggsave(muts_per_patient_bar, filename =  "../results/031_barplot_mutations_per_patient_cscape.png")
  
# 
  # DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
  #                                        into = "Patient_ID", sep = 12) %>%
  #   drop_na(Patient_ID) %>% 
  #   ggplot(aes(x= Patient_ID, fill = CScape_Mut_Class), alpha = 0.8) +
  #   geom_boxplot()
  

  
  
# Genes of interest ---------------------------
  # Gene with most driver muts
  top <- Oncogenic_mediators_mutation_summary %>%  
    filter(CScape_Driver >= 5) %>% pull(Hugo_Symbol)

  top_drivers <- DEG_Mutations_Annotations %>% filter(Hugo_Symbol %in% top) %>% 
    separate(col = Tumor_Sample_Barcode, into = "Patient_ID", sep = 12) %>% 
    ggplot(aes(y = Patient_ID, fill = CScape_Mut_Class), alpha = 0.8)+
    geom_bar()+
    facet_grid(cols = vars(Hugo_Symbol))
  ggsave(top_drivers, filename = '../results/031_barplot_top_drivers_in_patients.png')
  
  #SYNE1_ 
  DEG_Mutations_Annotations %>% filter(Hugo_Symbol == 'SYNE1') 
  DEG_Mutations_Annotations %>% filter(Hugo_Symbol == 'SYNE1') %>% 
    select(Start_Position, Reference, Mutant, Variant_Classification, CScape_Mut_Class,
           CScape_Coding_score, CScape_Noncoding_score)

  
  mymaf <- read.maf(maf = DEG_Mutations_Annotations)
  plotProtein(gene ='SYNE1', refSeqID = 'NM_182961')
  lollipopPlot(maf = mymaf, gene = 'SYNE1', AACol = 'HGVSp_Short', showMutationRate = TRUE, refSeqID = "NM_182961")
  lollipopPlot(maf = mymaf, gene = 'CD200', AACol = 'HGVSp_Short', showMutationRate = TRUE)

