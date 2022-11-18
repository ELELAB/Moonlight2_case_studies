# This script investigate all four subtypes including a survival analysis 
# This part is not used in Basal-like Case study ! 

# Loaf Libraries ----------------------------------------
library(tidyverse)
library(TCGAbiolinks)

# PART 1: Intristic subtypes -----------------------------

## Genes of interest ----------------
genelist <- c("ERBB2", "ESR1", "PGR") #Her2, Estrogen receptor, Progesterone repcetor


## Basal Tumor Samples ---------------
Basal_raw <- get(load("/data/raw_data/computational_data/TCGA_data_081220/Expression/cancer_subtypes/BRCA/BRCA.Basal/BRCA_BRCA.Basal_dataFilt_HUGO.rda"))
basal_tum <-TCGAquery_SampleTypes(barcode = colnames(Basal_raw),typesample = c("TP"))
col.num <- which(colnames(Basal_raw) %in% basal_tum)
Basal <- Basal_raw[,(col.num)] %>% as.data.frame() %>%  rownames_to_column(var = "Gene")

Basal_genes <- Basal %>% filter(Gene %in% genelist) %>% 
  pivot_longer(cols = -c("Gene"), names_to = "samples", values_to = "Expression") %>% 
  mutate(Type = "Tumor_Basal")

## LumA Tumor Samples ---------------
LumA_raw <- get(load("/data/raw_data/computational_data/TCGA_data_081220/Expression/cancer_subtypes/BRCA/BRCA.LumA/BRCA_BRCA.LumA_dataFilt_HUGO.rda"))
LumA_tum <-TCGAquery_SampleTypes(barcode = colnames(LumA_raw),typesample = c("TP"))
col.num <- which(colnames(LumA_raw) %in% LumA_tum)
LumA <- LumA_raw[,(col.num)] %>% as.data.frame() %>%  rownames_to_column(var = "Gene")

LumA_genes <- LumA %>% filter(Gene %in% genelist) %>% 
  pivot_longer(cols = -c("Gene"), names_to = "samples", values_to = "Expression") %>% 
  mutate(Type = "Tumor_LumA")

## LumB Tumor Samples ---------------
LumB_raw <- get(load("/data/raw_data/computational_data/TCGA_data_081220/Expression/cancer_subtypes/BRCA/BRCA.LumB/BRCA_BRCA.LumB_dataFilt_HUGO.rda"))
LumB_tum <-TCGAquery_SampleTypes(barcode = colnames(LumB_raw),typesample = c("TP"))
col.num <- which(colnames(LumB_raw) %in% LumB_tum)
LumB <- LumB_raw[,(col.num)] %>% as.data.frame() %>%  rownames_to_column(var = "Gene")

LumB_genes <- LumB %>% filter(Gene %in% genelist) %>% 
  pivot_longer(cols = -c("Gene"), names_to = "samples", values_to = "Expression") %>% 
  mutate(Type = "Tumor_LumB")

## Her2 Tumor Samples ---------------
Her2_raw <- get(load("/data/raw_data/computational_data/TCGA_data_081220/Expression/cancer_subtypes/BRCA/BRCA.Her2/BRCA_BRCA.Her2_dataFilt_HUGO.rda"))
Her2_tum <-TCGAquery_SampleTypes(barcode = colnames(Her2_raw),typesample = c("TP"))
col.num <- which(colnames(Her2_raw) %in% Her2_tum)
Her2 <- Her2_raw[,(col.num)] %>% as.data.frame() %>%  rownames_to_column(var = "Gene")

Her2_genes <- Her2 %>% filter(Gene %in% genelist) %>% 
  pivot_longer(cols = -c("Gene"), names_to = "samples", values_to = "Expression") %>% 
  mutate(Type = "Tumor_Her2")

## Normal Samples --------- 
Full_raw <- get(load("/data/raw_data/computational_data/TCGA_data_081220/Expression/cancer_types/BRCA/BRCA_dataFilt_HUGO.rda"))
Full_norm<- TCGAquery_SampleTypes(barcode=colnames(Full_raw), typesample = c("NT"))
col.num <- which(colnames(Full_raw) %in% Full_norm)
Full <- Full_raw[,(col.num)] %>% as.data.frame() %>%  rownames_to_column(var = "Gene")

Norm_genes <- Full %>% filter(Gene %in% genelist) %>% 
  pivot_longer(cols = -c("Gene"), names_to = "samples", values_to = "Expression") %>% 
  mutate(Type = "Normal")


## Join all data ------- 
all <-full_join(Basal_genes, Norm_genes) %>% 
full_join(LumA_genes) %>% 
  full_join(LumB_genes) %>% 
  full_join(Her2_genes)

## Plot ---------------------
ggplot(data =all, aes(y=Expression, x = Gene , fill = Type), alpha = 0.8)+
  geom_boxplot() +
  facet_wrap(vars(Gene), scales = "free")



# PART 2: ICH classification -----------------------------
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

all2 <- all %>% separate(col = samples, 
                         into = "Patient_ID", sep = 12, remove = FALSE) %>% 
  left_join(clinical, by = c("Patient_ID" = "submitter_id"))

#Barplot of count of patients in different stages colored to subtype
all2 %>% 
  ggplot(aes(x = ajcc_pathologic_stage, fill = Type), alphe = 0.8)+
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# PART 3: COX analysis -----------------------------------
library(survival)
library(survminer)

#Wrangle data --------------------------------------------

# Change variables to time and evet
data <- all2 %>%  select(Patient_ID, Type, days_to_last_follow_up, days_to_death,
                         vital_status, ajcc_pathologic_stage, age_at_index, age_at_diagnosis,
                         year_of_diagnosis) %>% 
  unique() %>% 
  mutate(Time = case_when(is.na(days_to_last_follow_up)~days_to_death,
                          !is.na(days_to_death)~days_to_death,
                          TRUE ~ days_to_last_follow_up ),
         Event = case_when(vital_status == "Alive"~0,
                           vital_status == "Dead"~1 )) %>% 
  filter(Type != 'Normal', 
         Time >= 0) %>% 
  mutate_if(is.character, as_factor) 

# check <-data %>% drop_na(age_at_diagnosis)
# check <-survfit(Surv( age_at_diagnosis, Time1, Event)~Type, data=check)
# plot(check, mark.time = TRUE, xscale = 365.25,xlab = "Years", ylab = "Survival", col = c(1,2,3,4))
# lLab <- gsub("x=","",names(check$strata)) 
# legend(
#   "topright",
#   legend=lLab,
#   col=1:4,
#   lty=1,
#   horiz=FALSE,
#   bty='n')

#Distribution of follow-up time.:
data %>% mutate(as_factor(Time)) %>% 
  ggplot(aes(x = Time, fill = as_factor(Event)))+
  geom_histogram(alpha = 0.8, position = "identity")

# Combine stages to only three stages in new data set 
data_stage <- data %>% mutate(Stage = case_when(ajcc_pathologic_stage  == 'Stage I'|
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
  mutate(age_group=cut(age_at_index,breaks=seq(10,120,by=30)))

#Distribution of Age: 
data_stage %>% 
  #group_by(age_group) %>% summarise(n=n())
  ggplot(aes(x = age_at_index), alpha = 0.8)+
  #geom_bar(aes(y = (..count..)/sum(..count..)))
  geom_density()+
  facet_wrap(vars(Stage))

#Are age distribution equal in the stages?
x <- data_stage %>% mutate_if(is.character, as_factor) %>% 
  select(Patient_ID, age_at_index, Stage) %>% filter(Stage != 'Unknown')
levels(x$Stage )
ggplot(x, aes(x = Stage, y = age_at_index, color = Stage), alpha = 0.8)+
  geom_violin()+
  geom_boxplot() 
aov(age_at_index~Stage, data = x)
summary(aov(age_at_index~Stage, data = x ))

#Are age distribution equal in the types?
x <- data_stage %>% mutate_if(is.character, as_factor) #%>% select(Patient_ID, age_at_index, Type)
levels(x$Type )
ggplot(x, aes(x = Type, y = age_at_index, color = Type), alpha = 0.8)+
  geom_boxplot()
#geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
#facet_wrap(vars(Stage))
aov(age_at_index~Type, data = x)
summary(aov(age_at_index~Type, data = x ))


# Survial analysis ---------------------------------------------
## STRATA on subtype and Age ----------------------------------------------

#viridis::viridis(4)
colors <- c("Tumor_Basal" = '#440154FF', "Tumor_Her2"= '#31688EFF', "Tumor_LumA" = '#35B779FF', "Tumor_LumB" = '#FDE725FF')

# Stratify only on subtype, plot 5 years and all years
sfit_t <- survfit(Surv(Time, Event)~Type, data = data_stage )
summary(sfit_t) 
summary(survfit(Surv(Time, Event)~Type, data = data_stage ), times = (365.25*5)) # What are the surv at five years?
# Basal = 79,5 % survives at least 5 years
# Luminal A = 86,5 %
# Luminal B = 76,8 %
# Her2 = 71,8 % 

ggsurvplot(sfit_t, color = "Type",# pval = TRUE, 
           palette = colors,
           xscale = "d_y", xlab = "Years", xlim = c(0,1800), break.time.by = 365) #, risk.table = TRUE) #, conf.int = TRUE)
ggsurvplot(sfit_t, color = "Type",# pval = TRUE, 
           palette = colors,
           xscale = "d_y", xlab = "Years", break.time.by = 365.25) #, risk.table = TRUE) #, conf.int = TRUE)
coxph(Surv(Time, Event)~Type, data = data_stage) %>% summary()


# Stratify on subtype and age at diagnosis, plot 5 years and all years 
sfit_ta <- survfit(Surv(Time, Event)~Type+age_group, data = data_stage )
summary(sfit_ta)
survdiff(Surv(Time, Event)~Type+age_group, data = data_stage)
ggsurvplot(sfit_ta, facet.by = c("age_group"), palette = colors,
           break.time.by = 730.5, xscale = "d_y", xlab = "Years")
ggsurvplot(sfit_ta, facet.by = c("age_group"), 
           palette = colors,
           xscale = "d_y", xlab = "Years", xlim = c(0,1800), break.time.by = 365) #, conf.int = TRUE)
coxph(Surv(Time, Event)~Type+age_group, data = data_stage) %>% summary()



# check: should we include stage? 
data_stage2 <- data_stage %>% filter(Stage != 'Unknown')
c1 <- coxph(Surv(Time, Event)~Type+age_at_index, data = data_stage)
c2 <- coxph(Surv(Time, Event)~Type+age_at_index+Stage, data = data_stage)
summary(c2)
anova(c1,c2, test = 'LRT') # There is a sig. difference with adding age 

# Check: can we drop age_group?
c1 <- coxph(Surv(Time, Event)~Type, data = data_stage)
c2 <- coxph(Surv(Time, Event)~Type+age_at_index, data = data_stage)
summary(c2)
anova(c1,c2, test = 'LRT') # There is a sig. difference with adding age 

# Check Assumptions: 
## Hazard ratio is proportional over time
cox.zph(c2)
par(mfrow =c(3,1))
plot(cox.zph(c2), 
     xlab = "Time in days",
     main = "Hazard Ratio Assumption Check")
abline(h=0, col=2)
# We want the avg to be as close to zero as possible, meaning no dependency on time

## Linearity of numerical Xs 
plot(predict(c2), residuals(c2, type = "martingale"),
     xlab = "fitted values",
     ylab = "Martingale res", main = "Residual plot", las = 1)
abline(h=0, col=2)
lines(smooth.spline(predict(c2),
                    residuals(c2, type ="martingale")),
      col = 'blue')


# Type, age and stage (without stage unknown)
#Subtype and age and stage
sfit_tas <- survfit(Surv(Time, Event)~Type+Stage, data = data_stage2 )
summary(sfit_tas)
survdiff(Surv(Time, Event)~Type+Stage, data = data_stage2)
ggsurvplot(sfit_tas, facet.by = c("Stage"), palette = colors,
           break.time.by = 730.5, xscale = "d_y", xlab = "Years")
ggsurvplot(sfit_tas, facet.by = c("Stage"), 
           palette = colors,
           xscale = "d_y", xlab = "Years", xlim = c(0,1800), break.time.by = 365) #, conf.int = TRUE)
coxph(Surv(Time, Event)~Type+age_group, data = data_stage) %>% summary()


