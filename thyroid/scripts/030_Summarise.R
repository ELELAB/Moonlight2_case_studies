# Libraries -------------
library(tidyverse)

# Data ------------------
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
DMA_Thyroid <- get(load("../data/DMA_Thyroid.rda"))

# Number of Patients ----------
# see 06_Metadata_analysis.R

# Numbers to Table---------------
drivers <- append(DMA_Thyroid$TSG, DMA_Thyroid$OCG)

drivers_transciption <- DEG_Mutations_Annotations %>% 
  filter(!is.na(Moonlight_Oncogenic_Mediator), 
         CScape_Mut_Class == 'Driver',
         Potential_Effect_on_Transcription == 1) %>% 
  pull(Hugo_Symbol) %>% unique()

## Without mutations
DEG_Mutations_Annotations %>% filter(!is.na(logFC),
                                     is.na(Start_Position)) %>% count()

Oncogenic_mediators_mutation_summary %>% 
  filter(Total_Mutations == 0) %>% group_by(Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

## Genes: 
Oncogenic_mediators_mutation_summary %>% 
  group_by(Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

Oncogenic_mediators_mutation_summary %>% 
  filter(Hugo_Symbol %in% drivers) %>% 
  group_by(Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

Oncogenic_mediators_mutation_summary %>% 
  filter(Hugo_Symbol %in% drivers_transciption) %>% 
  group_by(Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

## Mutations: 
DEG_Mutations_Annotations %>% 
  group_by(CScape_Mut_Class) %>% 
  summarise(n=n())

DEG_Mutations_Annotations %>% 
  filter(!is.na(Moonlight_Oncogenic_Mediator)) %>% 
  group_by(CScape_Mut_Class, Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

DEG_Mutations_Annotations %>% 
  filter(Hugo_Symbol %in% drivers) %>% 
  group_by(CScape_Mut_Class, Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

DEG_Mutations_Annotations %>% 
  filter(Hugo_Symbol %in% drivers_transciption) %>% 
  group_by(CScape_Mut_Class, Moonlight_Oncogenic_Mediator) %>% 
  summarise(n=n())

## Subset DEG_Mutations_Annotations to include only mutation information 
## of driver genes
driver_genes_mut <- DEG_Mutations_Annotations %>% 
  filter(Hugo_Symbol %in% drivers,
         CScape_Mut_Class == "Driver") %>% 
  dplyr::select(Chromosome, Start_Position, End_Position,
                Strand, Hugo_Symbol, logFC, adj.P.Val,
                contains("Moonlight"), contains("CScape"),
                contains("Effect"), contains("Annotation"),
                contains("NCG"), Variant_Classification,
                Variant_Type, contains("Allele"), contains("Barcode"),
                -c("all_effects"))
write_csv(driver_genes_mut, file = "../results/030_driver_genes_muts_annotations.csv")


# Boxplots Gene/s/Mutations 
# This plot is not meaningful  - something is wrong
DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
                                       into = "Patient_ID", sep = 12) %>%
  group_by(Patient_ID, Moonlight_Oncogenic_Mediator) %>% summarise(n= n()) %>% 
  ggplot(aes(x = n,  y= Moonlight_Oncogenic_Mediator))+
  geom_boxplot()


# Unclassified Mutations -------------------
# Chromosome X and Y 
DEG_Mutations_Annotations %>% filter(Chromosome == 'chrX', 
                                     !is.na(Moonlight_Oncogenic_Mediator)) %>% 
  group_by(Moonlight_Oncogenic_Mediator, Chromosome, Hugo_Symbol)  %>% 
  summarise(n =n())

#Unclassified
Oncogenic_mediators_mutation_summary %>% filter(is.na(CScape_Driver), 
                                                  is.na(CScape_Passenger),
                                                CScape_Unclassified >= 1) 
# Mutations in multiple patients ----------- 
DEG_Mutations_Annotations %>% 
  separate(col = Tumor_Sample_Barcode, 
           into = "Patient_ID", sep = 12) %>%
  select(Hugo_Symbol, Moonlight_Oncogenic_Mediator, Chromosome, Start_Position, End_Position, Variant_Classification, Variant_Type, CScape_Mut_Class) %>% 
  group_by(Hugo_Symbol, Moonlight_Oncogenic_Mediator, Chromosome, Start_Position, End_Position) %>% nest() %>% 
  mutate(n = map_dbl(data,nrow))  %>% 
  filter(n > 1) %>%  unnest(cols = c(data))
