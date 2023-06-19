library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(ggplot2)

# Load Data ---------------------------
clinical <- GDCquery_clinic(project = "TCGA-THCA", type = "clinical")
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
dataFilt_hugo <- get(load("../data/rawdata/THCA_dataFilt_HUGO.rda"))


# Extract patient IDs ---------------
IDs_mut <- DEG_Mutations_Annotations %>% separate(col = Tumor_Sample_Barcode, 
         into = "Patient_ID", sep = 12) %>% 
  pull(Patient_ID) %>% unique() %>% na.omit()

IDs_exp <- data.frame(id = colnames(dataFilt_hugo))%>%
  separate(col = id, into = "Patient_ID", sep = 12) %>% pull(Patient_ID) %>% unique()

IDs <- append(IDs_exp, IDs_mut) %>% unique()

# Which are unique to one data set: 
IDs_mut[!(IDs_mut %in% IDs_exp)] # one patient
IDs_exp[!(IDs_exp %in% IDs_mut)] # 20 patients


# Data Wrangle -------------------
clinical_basal <- clinical %>% filter(submitter_id %in% IDs)
colnames(clinical_basal)


# Plot metadata variables data ------------------------------
clinical_basal %>% 
  ggplot(aes(y = primary_diagnosis, fill = ajcc_pathologic_stage), alpha = 0.8) +
  geom_bar() +
  theme_bw()
  
clinical_basal %>% 
  ggplot(aes(y = primary_diagnosis, fill = race), alpha = 0.8) +
  geom_bar(position = "dodge") +
  theme_bw()

clinical_basal %>% 
  ggplot(aes(x = age_at_index, fill = race), alpha = 0.8) +
  geom_bar(position = 'dodge') +
  theme_bw()

clinical_basal %>% mutate(age_group=cut(age_at_index,breaks=seq(10,100,by=10))) %>% 
  ggplot(aes(x = age_group, fill = race), alpha = 0.8) +
  geom_bar() +
 #geom_bar(position = position_dodge2(preserve = "single")) +
  theme_bw()

clinical_basal %>% 
  ggplot(aes(y = age_at_index,x = ajcc_pathologic_stage, fill = ajcc_pathologic_stage), alpha = 0.8) +
  geom_boxplot() +
  theme_bw()

clinical_basal %>% 
  ggplot(aes(x = prior_malignancy)) +
  geom_bar()


# Combine Stages ------------------
data_stage <- clinical_basal %>% mutate(Stage = case_when(ajcc_pathologic_stage  == 'Stage I'|
                                                  ajcc_pathologic_stage  == 'Stage IA'| 
                                                  ajcc_pathologic_stage  == 'Stage IB'| 
                                                  ajcc_pathologic_stage  == 'Stage 1'|
                                                  ajcc_pathologic_stage  == 'Stage II'|
                                                  ajcc_pathologic_stage  == 'Stage IIA'~'Early Stage',
                                                ajcc_pathologic_stage  == 'Stage IIB'| 
                                                  ajcc_pathologic_stage  == 'Stage III'|
                                                  ajcc_pathologic_stage  == 'Stage IIIA'|
                                                  ajcc_pathologic_stage  == 'Stage IIIB'|
                                                  ajcc_pathologic_stage  == 'Stage IIIC'~'Localluy advanced', 
                                                ajcc_pathologic_stage  == 'Stage IV'~'Metastatic' ,
                                                ajcc_pathologic_stage  == 'Stage X' ~ 'Unknown',
                                                TRUE~'Unknown' )) %>% 
  mutate(age_group=cut(age_at_index,breaks=seq(10,100,by=10)))

#Distribution of Age across Stage: 
data_stage %>% 
  #group_by(age_group) %>% summarise(n=n())
  ggplot(aes(x = age_at_index), alpha = 0.8)+
  #geom_bar(aes(y = (..count..)/sum(..count..)))
  geom_density()+
  facet_wrap(vars(Stage))

#Are age distribution equal in the stages?
x <- data_stage %>% mutate_if(is.character, as_factor) %>% 
  select(age_at_index, Stage) %>% filter(Stage != 'Unknown')
levels(x$Stage )
ggplot(x, aes(x = Stage, y = age_at_index, color = Stage), alpha = 0.8)+
  geom_violin()+
  geom_boxplot() 
aov(age_at_index~Stage, data = x)
summary(aov(age_at_index~Stage, data = x ))


# Plotting with maftools -----------------

mymaf <- read.maf(maf = DEG_Mutations_Annotations)

oncoplot(mymaf, top = 20)
tcgaCompare(maf = mymaf, 
            cohortName = 'Entire Basal', 
            logscale = TRUE, 
            capture_size = 50)

rainfallPlot(maf = mymaf, 
             detectChangePoints = TRUE)
