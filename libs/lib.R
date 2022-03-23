library(tidyverse)
library(data.table)
# Read and process the emmaper results - RAW
read_emmaper_files <- function(file, pattern = pattern) {
  col_names <-
    c(
      "query_name",
      "seed_eggNOG_ortholog",
      "seed_ortholog_evalue",
      "seed_ortholog_score",
      "Predicted_taxonomic_group",
      "Predicted_protein_name",
      "Gene_Ontology_terms",
      "EC_number",
      "KEGG_ko",
      "KEGG_Pathway",
      "KEGG_Module",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "tax_scope",
      "eggNOG_OGs",
      "bestOG",
      "COG_Functional_Category",
      "eggNOG_free_text_description"
    )
  results <- fread(file,
        col.names = col_names,
        quote = "",
        sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(
      pattern = pattern,
      replacement = "",
      x = basename(file)
    )) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(
      eggNOG_OGs,
      into = c("og", "taxid"),
      sep = "@",
      remove = FALSE
    ) %>%
    left_join(taxmap, by = "taxid") %>%
    filter(Predicted_taxonomic_group == tax_name) %>%
    distinct() %>%
    left_join(annotations, by = "og") %>%
    mutate(bestOG = og)
  return(results)
}

# Read and process the emmaper results - LOWEST OG
read_emmaper_files_lower <- function(file, pattern = pattern) {
  col_names <-
    c(
      "query_name",
      "seed_eggNOG_ortholog",
      "seed_ortholog_evalue",
      "seed_ortholog_score",
      "Predicted_taxonomic_group",
      "Predicted_protein_name",
      "Gene_Ontology_terms",
      "EC_number",
      "KEGG_ko",
      "KEGG_Pathway",
      "KEGG_Module",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "tax_scope",
      "eggNOG_OGs",
      "bestOG",
      "COG_Functional_Category",
      "eggNOG_free_text_description"
    )
  results <- fread(file,
        col.names = col_names,
        quote = "",
        sep = "\t") %>%
    as_tibble() %>%
    mutate(sMAG = gsub(
      pattern = pattern,
      replacement = "",
      x = basename(file)
    )) %>%
    separate_rows(eggNOG_OGs, sep = ",") %>%
    separate(
      eggNOG_OGs,
      into = c("og", "taxid"),
      sep = "@",
      remove = FALSE
    ) %>%
    left_join(taxmap, by = "taxid") %>%
    group_by(sMAG, query_name) %>%
    #slice(which.min(taxid)) %>% View()
    filter(taxid == min(taxid)) %>%
    ungroup() %>%
    #filter(Predicted_taxonomic_group == tax_name) %>%
    distinct() %>%
    select(sMAG, query_name, og, taxid) %>%
    distinct() %>%
    mutate(bestOG = og) %>%
    inner_join(annotations %>% rename(taxid = level), by = c("og", "taxid")) %>%
    mutate(description = ifelse(description == "", "Not available", description))
  return(results)
}



filter_prev <- function(physeq, prev, data) {

  rare_taxa <- microbiome::rare_members(physeq, prevalence = prev)

  physeq_filt <- microbiome::remove_taxa(rare_taxa, physeq)

  coms <- taxa_names(physeq_filt)

  n_genes <- data %>%
    filter(comm_name %in% coms) %>%
    group_by(category) %>%
    count(name = "ngenes") %>%
    ungroup()

  n_genes <- n_genes %>%
    mutate(total_genes = n_genes$ngenes %>% sum())

  n_coms <- data %>%
    filter(comm_name %in% coms) %>%
    select(comm_name, category) %>%
    distinct() %>%
    group_by(category) %>%
    count(name = "ncoms") %>%
    ungroup()

  n_coms <- n_coms %>%
    mutate(total_coms = n_coms$ncoms %>% sum())

  res <- n_genes %>%
    inner_join(n_coms, by = "category") %>%
    select(category, ngenes, ncoms, total_genes, total_coms) %>%
    mutate(prev = prev)

  mag_summary <- data %>%
    group_by(sMAG) %>%
    count(name = "ngenes_all") %>%
    ungroup() %>%
    inner_join(data %>%
                 filter(comm_name %in% coms) %>%
                 group_by(sMAG) %>%
                 count(name = "ngenes_prev") %>%
                 ungroup(), by = "sMAG") %>%
    mutate(prop_genes_kept = ngenes_prev/ngenes_all,
           prev = prev)

  physeq_filt_df <- as(otu_table(physeq_filt), "matrix") %>%
    as.data.frame() %>%
    rownames_to_column("comm_names") %>%
    as_tibble()

  return(list(data = physeq_filt_df, cnts = res, mag_summary = mag_summary))
}


filter_prev_og <- function(physeq, prev, data) {

  rare_taxa <- microbiome::rare_members(physeq, prevalence = prev)

  physeq_filt <- microbiome::remove_taxa(rare_taxa, physeq)

  coms <- taxa_names(physeq_filt)

  n_genes <- data %>%
    filter(comm_name %in% coms) %>%
    group_by(category) %>%
    count(name = "ngenes") %>%
    ungroup()

  n_genes <- n_genes %>%
    mutate(total_genes = n_genes$ngenes %>% sum())

  n_coms <- data %>%
    filter(comm_name %in% coms) %>%
    select(comm_name, category) %>%
    distinct() %>%
    group_by(category) %>%
    count(name = "ncoms") %>%
    ungroup()

  n_coms <- n_coms %>%
    mutate(total_coms = n_coms$ncoms %>% sum())

  res <- n_genes %>%
    inner_join(n_coms) %>%
    select(category, ngenes, ncoms, total_genes, total_coms) %>%
    mutate(prev = prev)

  mag_summary <- data %>%
    group_by(sMAG) %>%
    count(name = "ngenes_all") %>%
    ungroup() %>%
    inner_join(data %>%
                 filter(comm_name %in% coms) %>%
                 group_by(sMAG) %>%
                 count(name = "ngenes_prev") %>%
                 ungroup()) %>%
    mutate(prop_genes_kept = ngenes_prev/ngenes_all,
           prev = prev)

  physeq_filt_df <- as(otu_table(physeq_filt), "matrix") %>%
    as.data.frame() %>%
    rownames_to_column("comm_names") %>%
    as_tibble()

  return(list(data = physeq_filt_df, cnts = res, mag_summary = mag_summary))
}


majority_vote <- function (x, seed = 12345) {
  set.seed(seed)
  whichMax <- function(x) {
    m <- seq_along(x)[x == max(x, na.rm = TRUE)]
    if (length(m) > 1)
      sample(m, size = 1)
    else m
  }
  x <- as.vector(x)
  tab <- table(x)
  m <- whichMax(tab)
  out <- list(table = tab, ind = m, majority = names(tab)[m])
  return(out)
}

apply_majority <- function(X){
  X <- setDT(X)
  DT.1 <- X[,majority:=majority_vote(bestOG)$majority, by="comm_name"]
  df <- DT.1 %>% as_tibble() %>% distinct()
}

get_majority <- function(X){

  list_genes <- X %>%                                        # Split into groups by gene-caller-id
    split(.$comm_name)

  maj_l <- future_map(list_genes,apply_majority, .progress = TRUE) # run majority_vote function
  maj_df <- plyr::ldply(maj_l, data.frame) %>%  # bind list rowwise and get distint votes for gene category
    select(comm_name,majority) %>%
    distinct() %>%
    as_tibble()
}
