# This script performs gene literature search of predicted driver genes

# GLS function definition. This was taken from commit 28205ac in the Moonlight2R
# repository. This was done because we want to keep the old Moonlight2R
# version, that we actually used to generate the data, in this workflow
# and environment
GLS <- function(genes,
                query_string = "AND cancer AND driver",
		max_records = 20) {

  # Initialize empty tibble to store results
  pubmed_mining <- tibble()

  # For each gene x in input, search PubMed based on specified
  # query
  literature_search <- map(genes, function(x) {

    pubmed_query <- paste(x, query_string)

    # Search and retrieve results from PubMed
    gene_pubmed <- get_pubmed_ids(pubmed_query)

    # Retrieve number of publications
    count_pubmed <- gene_pubmed$Count %>%
      as.numeric()

    # If query matches any pubmed records
    if (count_pubmed > 0) {

      # Fetch data of PubMed records searched via above query
      top_results <- fetch_pubmed_data(gene_pubmed,
                                       retstart = 0,
                                       retmax = max_records)

      # Extract information from PubMed records into a table
      record_info <- table_articles_byAuth(top_results,
                                           included_authors = "first",
                                           max_chars = -1,
                                           getKeywords = TRUE)

      # Select only PubMed id, doi, title, abstract, year, and keywords of
      # PubMed records
      record_info_wrangled <- record_info %>%
        as_tibble() %>%
        dplyr::select(c(pmid, doi, title, abstract, year, keywords)) %>%
        mutate(gene = x,
               pubmed_count = count_pubmed) %>%
        dplyr::relocate(gene,
                        .after = pmid)

      # Bind table to table containing results from previous gene(s)
      pubmed_mining <- pubmed_mining %>%
        bind_rows(record_info_wrangled)

      # If no records of query is found in PubMed
    } else {

      # Create tibble of one row of gene that did not have any PubMed results
      no_results_tbl <- tibble(pmid = NA,
                               gene = x,
                               doi = NA,
                               title = NA,
                               abstract = NA,
                               year = NA,
                               keywords = NA,
                               pubmed_count = count_pubmed)

      # Bind tibble of gene without PubMed information to table containing
      # results of previous gene(s)
      pubmed_mining <- pubmed_mining %>%
        bind_rows(no_results_tbl)

    }

  }) %>%
    bind_rows()

  return(literature_search)

}

# Load data
DMA_Basal <- get(load("../data/DMA_Basal.rda"))

# Wrangle data
drivers <- append(DMA_Basal$TSG, DMA_Basal$OCG)

# Perform literature search of driver genes
GLS_data <- GLS(genes = drivers)

# Save results
write_csv(GLS_data, file = "../results/08_driver_literature.csv")

