# This script compare known lung cancer driver genes with moonlight findings.

# SET UP -------------------------------
# Library ----------------
library(tidyverse)
library(Moonlight2R)

# Data -------------------
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
GUST <- read.csv("../data/rawdata/GUST_LUAD.csv")
data(NCG)

# Analysis -------------------------
Drivers <- Oncogenic_mediators_mutation_summary %>% filter( CScape_Driver >= 1)

## Extract NCG lung cancer candidates ----------------------
NCG %>% filter(str_detect(NCG_cancer_type,'lung')) %>% count()

## Oncogenic mediators which are NCG lung cancer candidates
breast_cand <- Oncogenic_mediators_mutation_summary %>% 
  filter(str_detect(NCG_cancer_type,'lung')) %>% pull(Hugo_Symbol)

## Driver mutations in above genes 
Drivers %>% 
  filter(Hugo_Symbol %in% breast_cand)

# GUST vs Moonlight Breast cancer ----------------------------
full_join(Drivers, GUST, by = c("Hugo_Symbol" = "Symbol")) %>% 
  group_by(Moonlight_Oncogenic_Mediator, GUST_pred) %>% 
  summarise(n = n())

moonlight_gust_tab <- left_join(Drivers, GUST, by = c("Hugo_Symbol" = "Symbol")) %>% 
  filter(!is.na(GUST_pred)) %>% 
  select(Hugo_Symbol, Moonlight_Oncogenic_Mediator, GUST_pred, NCG_driver, NCG_cancer_type) 

write.csv(moonlight_gust_tab, file = "../results/07_overlapping_genes_moonlight_ncg_gust.csv")

## CScape-somatic threshold score -----------------------------
# Should we adjust driver threshold based on NGC and CScape-scores, as 
# suggested by the paper. 
pivot <- DEG_Mutations_Annotations %>% 
  pivot_longer(cols = c(CScape_Coding_score, CScape_Noncoding_score), 
               names_to = "Type", values_to =  'CScape_score') # %>%
specific <- pivot %>%  filter( Hugo_Symbol %in% breast_cand)

ggplot(data = pivot, aes(x = CScape_score, y = Type))+
  geom_boxplot() +
  geom_jitter() +
  geom_jitter(data = specific, aes(x = CScape_score, y = Type, color = CScape_Mut_Class), size = 3)
#Conclusion: no, if they had lumped together with high scores then we should. 

