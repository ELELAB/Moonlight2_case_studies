# An enrichment analysis of moonlights predictions of OCG and TSG

# Libraries -------------------------
library(tidyverse)
library(enrichR)
library(gridExtra)
source("99_functions.R")

#Load data -------------------------
Oncogenic_mediators_mutation_summary <- get(load("../results/Oncogenic_mediators_mutation_summary.rda"))
DEG_Mutations_Annotations <- get(load("../results/DEG_Mutations_Annotations.rda"))
DMA_Lung <- get(load("../data/DMA_Lung.rda"))

# Kinases  & Transcription factors --------------------------------
kinhub <- read_tsv("../data/rawdata/kinhub.tsv")




# Mutation Comparison ----------------------------------------------
DEG_Mutations_Annotations %>% filter(!is.na(Moonlight_Oncogenic_Mediator), 
                                     CScape_Mut_Class != 'No_mutations') %>% 
  ggplot(aes(x = Variant_Classification, fill = CScape_Mut_Class),alpah = 0.8)+
  geom_bar()+
  facet_wrap(facets = vars(Moonlight_Oncogenic_Mediator)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#all mutations in ocg and tsg in percentage
putative_var_class<- DEG_Mutations_Annotations %>% filter(!is.na(Moonlight_Oncogenic_Mediator),
                                                          CScape_Mut_Class != 'No_mutations') %>% 
  ggplot(aes(x = Variant_Classification) , alpah = 0.8)+
  geom_bar(aes(y = ..prop.., group = 1, fill = 1), stat = 'count')+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(facets = vars(Moonlight_Oncogenic_Mediator)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(putative_var_class, filename = "../results/05_variant_class_percent_putative_TSG_OCG.png")

#same as above but only driver mutations 
driver_var_class <-DEG_Mutations_Annotations %>% filter(!is.na(Moonlight_Oncogenic_Mediator),
                                                        #CScape_Mut_Class != 'No_mutations'.
                                                        CScape_Mut_Class == 'Driver') %>% 
  ggplot(aes(x = Variant_Classification) , alpah = 0.8)+
  geom_bar(aes(y = ..prop.., group = 1, fill = 1), stat = 'count')+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(facets = vars(Moonlight_Oncogenic_Mediator)) +
  theme_bw()+
  labs(y = 'Percent',  x = 'Variant Classification')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(driver_var_class, filename = "../results/05_variant_class_percent_driver_TSG_OCG.png")
ggsave(driver_var_class, filename = "../results/05_variant_class_percent_driver_TSG_OCG.pdf",
       width = 15, height = 12, units = "cm" )





# Enrichment Analysis ----------------------------------------------
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", 
         "GO_Biological_Process_2021", "KEGG_2021_Human")

Drivers <- Oncogenic_mediators_mutation_summary %>% filter( CScape_Driver >= 1)
Moonlight_TSG <- DMA_Lung$TSG 
Moonlight_OCG <- DMA_Lung$OCG 

oncogenes_enrich <- enrichr(genes = Moonlight_OCG, databases = dbs)
tumor_enrich <- enrichr(genes =Moonlight_TSG, databases = dbs)

p2 <- goplot(data = oncogenes_enrich$GO_Molecular_Function_2021, title = "Oncogenes" , top =10)
p3 <- goplot(data = tumor_enrich$GO_Molecular_Function_2021, title = "Tumor suppressors", top =10)
grid.arrange(p2, p3,ncol=2)

p2 <- goplot(data = oncogenes_enrich$GO_Biological_Process_2021, title = "Oncogenes" , top =10)
p3 <- goplot(data = tumor_enrich$GO_Biological_Process_2021, title = "Tumor suppressors", top =10)
grid.arrange(p2, p3,ncol=2)

p2 <- goplot(data = oncogenes_enrich$GO_Cellular_Component_2021, title = "Oncogenes" , top =10)
p3 <- goplot(data = tumor_enrich$GO_Cellular_Component_2021, title = "Tumor suppressors", top =10)
grid.arrange(p2, p3,ncol=2)

p2 <- goplot(data = oncogenes_enrich$KEGG_2021_Human, title = "Oncogenes" , top =10)
p3 <- goplot(data = tumor_enrich$KEGG_2021_Human, title = "Tumor suppressors", top =10)
grid.arrange(p2, p3,ncol=2)

### All in one plot -----------------------------------------
onco_complete <- rbind(oncogenes_enrich$GO_Molecular_Function_2021, oncogenes_enrich$GO_Biological_Process_2021,
                       oncogenes_enrich$KEGG_2021_Human) %>% filter(Adjusted.P.value < 0.05) 
tumor_complete <- rbind( tumor_enrich$GO_Molecular_Function_2021, tumor_enrich$GO_Biological_Process_2021,
                         tumor_enrich$KEGG_2021_Human) %>% filter(Adjusted.P.value < 0.05)
p2 <- goplot(data = onco_complete, title = "Oncogenes" , top =15)
p3 <- goplot(data = tumor_complete, title = "Tumor suppressors", top =15)
go_plot <- grid.arrange(p2, p3,ncol=2 )
ggsave(go_plot, filename = "../results/05_GO_KEGG_top_adjpval_OCG_TSG.tiff",
       width = 19, height = 12, units = "cm", dpi = 300, compression = "lzw")
ggsave(go_plot, filename = "../results/05_GO_KEGG_top_adjpval_OCG_TSG.pdf",
       width = 19, height = 12, units = "cm")


#Share any significant terms?
inner_join(onco_complete, tumor_complete, by = c("Term"))



### Top 3 from each category ------------------------------------------
onco_complete <- rbind(slice_min(oncogenes_enrich$GO_Molecular_Function_2021,order_by = Adjusted.P.value, n =3), 
                       slice_min(oncogenes_enrich$GO_Biological_Process_2021,order_by = Adjusted.P.value, n =3),
                       slice_min(oncogenes_enrich$KEGG_2021_Human,order_by = Adjusted.P.value, n =3))
tumor_complete <- rbind(slice_min(tumor_enrich$GO_Molecular_Function_2021,order_by = Adjusted.P.value, n =3),
                        slice_min(tumor_enrich$GO_Biological_Process_2021,order_by = Adjusted.P.value, n =3),
                        slice_min(tumor_enrich$KEGG_2021_Human,order_by = Adjusted.P.value, n =3))

p2 <- goplot(data = onco_complete, title = "Oncogenes" , top =10)
p3 <- goplot(data = tumor_complete, title = "Tumor suppressors", top =10)
grid.arrange(p2, p3,ncol=2)


