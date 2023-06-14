# Functions 

# Mona's Upset plot functions ------------------------
library(UpSetR)
library(ComplexHeatmap)

### This function, upset_plot, creates an UpSet plot of intersections between
### different sets. This function takes five inputs:
### file_name:      File name of UpSet plot to be saved as a pdf
### dir_output:     Directory where UpSet plot should be saved
### sets_list:      A list containing those sets that are to be compared 
### names_sets:     A character vector containing names of each element in the list, sets_list
### title_plot:     The title to be included in the UpSet plot
### The function returns an UpSet plot. 
upset_plot <- function(file_name,
                       dir_output,
                       sets_list, 
                       names_sets, 
                       title_plot, viridis_color = "turbo") {
  
  # Save UpSet plot as pdf
  # tiff(file = paste0(dir_output, "/", file_name),
  #      width = 18, height = 8, units = "in", pointsize = 12,
  # compression ="lzw", res = 300)
  pdf(file = paste0(dir_output,
                    "/",
                    file_name), compress = FALSE,
      width = 18,
      height = 8)
  
  # Add names to elements in list of sets
  names(sets_list) <- names_sets
  
  # Make combination matrix of sets to be used in UpSet plot
  comb_mat_sets <- make_comb_mat(sets_list)
  
  # Generate color palette for plot using viridis package
  n_col <- max(comb_degree(comb_mat_sets))
  palette_col <- viridis_pal(option = viridis_color)(n_col)
  
  # Create UpSet plot
  upset_p <- UpSet(comb_mat_sets, 
                   set_order = names(sets_list),
                   pt_size = unit(5, 
                                  "mm"), 
                   lwd = 3, 
                   height = unit(4, 
                                 "cm"),
                   comb_col = palette_col[comb_degree(comb_mat_sets)],
                   top_annotation = upset_top_annotation(comb_mat_sets, 
                                                         height = unit(12, 
                                                                       "cm"),
                                                         ylim = c(0, 
                                                                  max(comb_size(comb_mat_sets))),
                                                         bar_width = 0.7, 
                                                         axis_param = list(side = "left", 
                                                                           at = seq(from = 0,
                                                                                    to = max(comb_size(comb_mat_sets)),
                                                                                    by = 500)),
                                                         annotation_name_side = "left", 
                                                         annotation_name_gp = gpar(cex = 1), 
                                                         annotation_name_offset = unit(1.5,
                                                                                       "cm")),
                   right_annotation = upset_right_annotation(comb_mat_sets, 
                                                             width = unit(3, 
                                                                          "cm"), 
                                                             gp = gpar(fill = "darkseagreen"),
                                                             axis_param = list(at = seq(from = 0,
                                                                                        to = max(set_size(comb_mat_sets)),
                                                                                        by = 2000)), 
                                                             annotation_name_offset = unit(1.5, 
                                                                                           "cm")),
                   row_names_gp = gpar(fontsize = 12))
  
  # Add number of elements in each set on top of bars in plot
  draw_upset <- draw(upset_p)
  col_ord <- column_order(draw_upset)
  c_s <- comb_size(comb_mat_sets)
  decorate_annotation("intersection_size", {
    grid.text(c_s[col_ord], 
              x = seq(c_s), 
              y = unit(c_s[col_ord], 
                       "native") + 
                unit(2, "pt"), 
              gp = gpar(fontsize = 12, 
                        fontface = "bold"),
              just = "bottom",
              default.units = "native")
  })
  
  # Add title to plot
  grid.text(label = title_plot, 
            x = unit(20, 
                     "cm"), 
            y = unit(18, 
                     "cm"), 
            gp = gpar(fontsize = 18),
            just = "centre")
  
  # End with dev.off() to save plot
  dev.off()
  
}


# Go Enrichemnt Point plot --------------------------------------------
library(enrichR)

goplot <- function(data, title = c(""), top = 15){
  myplot <- data %>% separate(col = Overlap, into = c('Count','Total'), 
                              sep ="/", remove =FALSE) %>% 
    mutate(Count  = as.numeric(Count),
           Gene.Ratio = round(Count/as.numeric(Total), 2)) %>% 
    slice_min(order_by = Adjusted.P.value,n = top, with_ties = FALSE) %>% 
    ggplot(aes(x = Gene.Ratio, y = reorder(Term, -Adjusted.P.value), color = Adjusted.P.value , size = Count )) +
    geom_point() +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width =30))+
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + 
    #theme(text = element_text(size = 16))+
    theme(
      text = element_text(size = 8))+#, family="Arial"))+
      #legend.key.size = unit(0.4, "cm")) +
      #axis.title.x = element_text(size= 4),#, family="Arial"),
      #axis.title.y = element_text(size= 4))+#, family="Arial"))+
    ylab("") + 
    xlab("Gene Ratio") + 
    ggtitle(title) +
    guides(
      color = guide_colorbar(order = 1),
      fill = guide_legend(order = 1))
  return(myplot)
}


# plotMoonlight
plotMoonlight2 <- function(DEG_Mutations_Annotations, 
                          Oncogenic_mediators_mutation_summary,
                          dataURA,
                          gene_type = "drivers",
                          n = 50, 
                          genelist = c(),
                          additionalFilename = ""){
  
  # The differentially expressed genes, that are annotated as TSG/OCG
  DEGs <- DEG_Mutations_Annotations %>% 
    select(Hugo_Symbol, Moonlight_gene_z_score, logFC) %>% 
    unique() %>% 
    drop_na(Moonlight_gene_z_score)
  
  # restructure URA to tibble
  ura <- as_tibble(dataURA, rownames = NA) %>% 
    rownames_to_column(var = "Genes")
  
  ura_wrangled <- ura %>%  
    pivot_longer(cols = !c('Genes'), 
                 names_to = 'Biological_Process',
                 values_to = 'Moonlight_score') %>% 
    right_join(Oncogenic_mediators_mutation_summary, 
               by = c("Genes" = "Hugo_Symbol")) %>% 
    right_join(DEGs, by = c("Genes" = "Hugo_Symbol")) %>%
    replace_na(list(CScape_Driver = 0, 
                    CScape_Passenger = 0, 
                    CScape_Unclassified = 0)) #%>% 
  
  # Type of plot:
  if (length(genelist) > 0 ){
    ura_wrangled <- ura_wrangled %>% 
      filter(Genes %in% genelist)
    n <- ura_wrangled %>% group_by(Moonlight_Oncogenic_Mediator) %>% summarise() %>% count() %>% pull()
    if (n <= 1){
      stop("The genelist must contain at least one OCG and one TSG.")
    }
    
  } else if (gene_type == "mediators"){
    ura_wrangled <- ura_wrangled %>% 
      slice_max(Total_Mutations, n = n, with_ties = FALSE)
    
  } else{
    ura_wrangled <- ura_wrangled %>% 
      slice_max(CScape_Driver, n = n, with_ties = FALSE)
    
  }
  
  # Variable for color scaling in legend
  max_driver <- ura_wrangled %>% arrange(desc(CScape_Driver)) %>% 
    select(CScape_Driver) %>%  head(1) %>% pull
  
  # Plot Heatmap
  bp_heatmap <- heatmap(ura_wrangled,
                        .row = Biological_Process,
                        .column = Genes,
                        .value = Moonlight_score,
                        scale = "none",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "complete",
                        cluster_rows = FALSE) %>%
    add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>%
    add_tile(logFC, palette = c("chartreuse4","firebrick3")) %>%
    #add_tile(CScape_Driver, palette = c(lower_col, "steelblue4")) %>% 
    #Use below when tidyHeatmap is updated
    #add_tile(CScape_Driver, palette = colorRamp2(c(0,max_driver), c("white", "dodgerblue3"))) %>%
    add_bar(CScape_Driver) %>% # Because it can happen that CScape_Driver can only have 1 distinct value  
    add_bar(Total_Mutations) 
  
  save_pdf(bp_heatmap, height = 15, width = 35, units = "cm",
           filename = paste(additionalFilename,"moonlight_heatmap.pdf", sep =""))
}
