# This scriptcontisns the investigation of MAF quality 
# the purpose were to determine thresholds of VAF, tumor and normal depth, and 
# variant alle count. 

# Libraries -----------------------
library(tidyverse)
library(gridExtra)

# Load data -----------------------
maf <- read.csv("../data/rawdata/mutations.csv")  %>% 
  mutate(VAF = t_alt_count/t_depth,
                Patient_ID = as.factor(X ))


# Investigation -------------------
# General density plots of scores ------------------------------
maf %>% filter(t_depth < 600) %>% 
  ggplot(aes(x = t_depth))+
  geom_density()

maf %>% #filter(n_depth < 600) %>% 
  ggplot(aes(x = n_depth))+
  geom_density()

maf %>% #filter(t_alt_count < 200) %>% 
  ggplot(aes(x = t_alt_count))+
  geom_density()


maf  %>% filter(VAF < 100) %>% 
  ggplot(aes(x = VAF)) +
  geom_density()

maf %>% 
  ggplot(aes(x = Patient_ID)) +
  geom_bar()


# ALL in one plot: 
# maf %>% mutate(VAF = t_alt_count/t_depth) %>% 
#   filter(VAF < 0.25) %>% arrange(VAF) %>% 
#   muatete(ID = row_number(-Score))
#   pivot_longer(cols = c(t_alt_count, t_depth, VAF), 
#                 names_to = 'Type',
#               values_to = 'Number') %>% 
#   ggplot(aes(x = fct_reorder2(.f =ID, VAF), y = Number, color = Type)) +
#   geom_point() +
#   facet_wrap(vars(Type), scale = 'free', nrow = 3)

#split tree way sorted by VAF -----------------------------------------
p1 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, VAF))) +
  geom_point(aes(y = VAF))

p2 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, VAF))) +
  geom_point(aes(y=t_depth))
p3 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, VAF))) +
  geom_point(aes(y=t_alt_count))
gridExtra::grid.arrange(p1,p2,p3, nrow = 3)

#split tree way sorted by t_alt_count -------------------------------------
pp1 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                      ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y = VAF))

pp2 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                      ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_depth))
pp3 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                      ID = as.factor(...1 )) %>% 
  filter(VAF < 0.25,
         t_alt_count >= 3,
         t_depth >= 10) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_alt_count))
gridExtra::grid.arrange(pp1,pp2,pp3, nrow = 3)


# Count mutations --------------------------------------------
#Mutations with low VAF and Count 
maf %>% filter(VAF > 0.05, #     %>% count()#,
               t_alt_count >= 3) %>% count()
#t_depth >= 10) 
maf %>% filter(VAF < 0.05, t_depth < 30 ) %>% count()

maf %>% filter(VAF < 0.05 | t_depth < 50 | t_alt_count <= 3) %>% count()

# Check 
p1 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.05,
         t_alt_count >= 8) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y = VAF))

p2 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.05,
         t_alt_count >= 8) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_depth))

p3 <- maf %>% mutate(VAF = t_alt_count/t_depth,
                     ID = as.factor(...1 )) %>% 
  filter(VAF < 0.05,
         t_alt_count >= 8) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_alt_count))

gridExtra::grid.arrange(p1,p2,p3, nrow = 3)






# Sort Patients into groups ----------------------------------------------- 
high_patients <- maf %>% separate(col = Tumor_Sample_Barcode,
                                  into = "Patient_ID", sep = 12) %>%
  group_by(Patient_ID) %>%  summarise(n = n()) %>%
  slice_max(n, n = 7) %>% filter(!is.na(Patient_ID)) %>% pull(Patient_ID)

low_patients <- maf %>% separate(col = Tumor_Sample_Barcode,
                                 into = "Patient_ID", sep = 12) %>%
  filter(!(Patient_ID %in% high_patients)) %>%  pull(Patient_ID) %>% unique()

maf_p <- maf %>%  separate(col = Tumor_Sample_Barcode,
                           into = "Patient_ID", sep = 12) %>% 
  mutate(patient_group = case_when(Patient_ID %in% high_patients ~ 'High',
                                   Patient_ID %in% low_patients ~'Low'))


## Investigate difference between patient groups ---------------- 
maf_p %>% filter(VAF > 0.05 | t_depth > 50 | t_alt_count >= 5)  %>% 
  ggplot(aes(x = VAF, color = patient_group), alpha = 0.8) +
  geom_density()

maf_p %>% filter(t_depth < 400) %>% 
  ggplot(aes(x = t_depth, color = patient_group), alpha = 0.8)+
  geom_density()

maf_p %>% filter(n_depth < 400) %>% 
  ggplot(aes(x = n_depth, color = patient_group), alpha = 0.8)+
  geom_density()

maf_p %>% filter(t_alt_count < 300) %>% 
  ggplot(aes(x = t_alt_count,color = patient_group), alpha = 0.8)+
  geom_density()


maf_p %>% #filter(VAF < 0.05 | t_depth < 30 | t_alt_count <= 5) %>% 
  ggplot(aes(x = patient_group, y = VAF, color = patient_group), alpha = 0.9)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))

maf_p %>% filter(VAF > 0.05 | t_depth > 30 | t_alt_count >= 5)  %>% 
  ggplot(aes(x = patient_group, y = VAF, color = patient_group), alpha = 0.9)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))

x1 <- maf_p %>% filter(VAF < 0.2 | t_depth < 10 | t_alt_count <= 3, 
                       patient_group == 'High') %>%
  pull(VAF)
x2 <- maf_p %>% filter(VAF < 0.2 | t_depth < 10 | t_alt_count <= 3,
                       patient_group == 'Low') %>% pull(VAF)
t.test(x1,x2, alternative = 'two.sided')


plot1 <- maf_p %>%  
  group_by(Patient_ID, patient_group) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Patient_ID, y = count, fill = patient_group)) +
  geom_col() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_text(aes(label=ifelse(patient_group=="High", 
                             Patient_ID, 
                             "")))+ 
  scale_y_continuous(limits = c(0,1900))

plot2 <- maf_p %>%  filter(VAF > 0.05 & n_depth > 10 & t_depth > 30, t_alt_count >= 5)   %>% 
  group_by(Patient_ID, patient_group) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Patient_ID, y = count, fill = patient_group)) +
  geom_col() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_text(aes(label=ifelse(patient_group=="High", 
                             Patient_ID, 
                             ""))) + 
  scale_y_continuous(limits = c(0,1900))

plot1_2<-gridExtra::grid.arrange(plot1,plot2, nrow = 2)
ggsave(plot1_2, filename = "../results/98_mutations_per_patient.png")
ggsave(plot1_2, filename = "../results/98_mutations_per_patient.pdf",
       width = 18, height = 12, units = "cm")


maf_p %>% filter(VAF > 0.05 & n_depth > 10 & t_depth > 30, t_alt_count >= 5)  %>% 
  ggplot() +
  geom_density(aes(x = VAF, color = Patient_ID),alpha = 0.8)+
  theme(legend.position = "none")



# FOUR PATIENTS 
# These patients were still outlieres in terms of mutations after filtering

four_patients<- c("TCGA-AO-A128","TCGA-BH-A18G","TCGA-D8-A1XK", "TCGA-D8-A1XQ")
mini_maf <- maf_p %>% filter(Patient_ID %in% four_patients) %>% 
  filter(VAF > 0.05 & n_depth > 10 & t_depth > 30, t_alt_count >= 5) 

mini_maf %>% ggplot(aes(x = VAF, color = Patient_ID), alpha = 0.8) +
  geom_density()

p1 <- mini_maf %>% mutate(VAF = t_alt_count/t_depth,
                          ID = as.factor(...1 )) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y = VAF))

p2 <- mini_maf %>% mutate(VAF = t_alt_count/t_depth,
                          ID = as.factor(...1 )) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_depth))

p3 <- mini_maf %>% mutate(VAF = t_alt_count/t_depth,
                          ID = as.factor(...1 )) %>% 
  ggplot(aes(x = fct_reorder(.f =ID, t_alt_count))) +
  geom_point(aes(y=t_alt_count))

gridExtra::grid.arrange(p1,p2,p3, nrow = 3)


# Conclusion ---------------------
# VAF > 0.05 
# t_depth > 30
# n_depth > 10
# t_alt_count >= 3 
