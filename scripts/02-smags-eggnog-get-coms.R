library(tidyverse)
library(data.table)
library(phyloseq)
library(microbiome)
library(furrr)
library(doFuture)
source("libs/lib.R")
options(future.globals.maxSize=10485760000)
registerDoFuture()
# You might want to change it
ncores <- 2
plan(future::multisession, workers = ncores)

# Table cleaning ----------------------------------------------------------

# read agnostos results
all_cls <- fread("data/agnostos/euk_all_genes_clusters_communities.tsv.gz", header = T)
all_cls <- all_cls %>%
  separate(original_gene_name, sep = " ", into = "query_name", extra = "drop", remove = FALSE) %>%
  #rename(sMAG = smag) %>%
  mutate(sMAG = gsub("mRNA.", "", smag)) %>%
  mutate(sMAG = gsub("_CNS_annotations-pep", "", sMAG),
         query_name = gsub(".+_CNS_annotations-pep_", "", query_name)) %>%
  as_tibble()

snglts <- read_tsv("data/agnostos/euk_singletons.tsv.gz")

# read the files produced in step 1
files <- list.files("results/lowestOG",pattern="*-filt.tsv.gz", full.names=TRUE)
annot <- data.table::rbindlist(future_map(files, fread, .progress = TRUE))
annot <- annot %>%
  mutate(sMAG = gsub("_CNS_annotations-pep", "", sMAG)) %>%
  as_tibble()



# Some sanity checks
# 4485410
annot %>% mutate(comb = paste0(query_name,"##",sMAG)) %>% .$comb %>% uniqueN()

# We will remove anything classified as MGE
to_remove <- annot %>%
  select(query_name, bestOG, description, sMAG) %>%
  filter(description == "transposition" | description == "Transposase IS4" | description == "Ribonuclease H") %>% as_tibble()

# 96395
to_remove %>% mutate(comb = paste0(query_name,"##",sMAG)) %>% .$comb %>% uniqueN()

annot1 <- annot %>%
  mutate(comb = paste0(query_name,"##",sMAG)) %>%
  filter(!(comb %in% (to_remove %>% mutate(comb = paste0(query_name,"##",sMAG)) %>% .$comb))) %>%
  select(-comb)

annot1 %>% mutate(comb = paste0(query_name,"##",sMAG)) %>% .$comb %>% uniqueN()


annot %>% left_join(all_cls) %>% as_tibble() %>% filter(is.na(cl_name)) %>% select(sMAG, query_name) %>% distinct()

annot %>% filter(grepl("Parent=Gene1", query_name)) %>% filter(sMAG == "TOSAG00-1") %>% select(sMAG, query_name)
all_cls %>% filter(grepl("Parent=Gene1", query_name)) %>% filter(sMAG == "TOSAG00-1") %>% select(sMAG, query_name)

annot %>% filter(grepl("TOSAG00-1", sMAG)) %>% select(sMAG) %>% distinct()
all_cls %>% filter(grepl("TOSAG00-1", sMAG)) %>% select(sMAG) %>% distinct()


# Total number of genes: 10,207,435
# How many genes where annotate by eggnog: 4,485,410
# Genes removed: 96,395
# In clusters:
# Discarded: 575053
#   with eggNOG: 453,283
# Singletons: 4,264,489
#   with eggNOG: 881,227
# Good: 5,367,893
#   with eggNOG: 3,054,505
#   without eggNOG: 2,313,388


all_cls %>% group_by(category) %>% count()

all_cls_annot <- all_cls %>% left_join(annot) %>% as_tibble()


disc <- all_cls_annot %>% filter(category == "DISC")
# 575053
# disc$original_gene_name %>% uniqueN()
disc_annot <- disc %>% filter(!is.na(bestOG))
# 469024
# disc_annot$original_gene_name %>% uniqueN()
sngl <- all_cls_annot %>% filter(category == "SINGL")
# 4264489
# sngl$original_gene_name %>% uniqueN()
sngl_annot <- sngl %>% filter(!is.na(bestOG))
# 892677
# sngl_annot$original_gene_name %>% uniqueN()

cls <- all_cls_annot %>% filter(category != "DISC", category != "SINGL")
# 5367893
# cls$original_gene_name %>% uniqueN()
cls_annot <- cls %>% filter(!is.na(bestOG))
# 3123709
# cls_annot$original_gene_name %>% uniqueN()
cls_noannot <- cls %>% filter(is.na(bestOG))
# 2244184
# cls_noannot$original_gene_name %>% uniqueN()

to_remove <- cls_annot %>%
  filter(description == "transposition" | description == "Transposase IS4" | description == "Ribonuclease H") %>% as_tibble()


# We will create tables that can be imported into anvi'o

comms <- all_cls %>%
  filter(!is.na(comm_name)) %>%
  filter(!(comm_name %in% (to_remove$comm_name %>% unique())))

comms_n <- comms %>%
  select(sMAG, comm_name) %>%
  group_by(sMAG, comm_name) %>%
  count(sort = TRUE) %>%
  ungroup()

# Convert to a data.frame
comms_df <- comms_n %>%
  spread(key = sMAG, value = n, fill = 0) %>%
  column_to_rownames(var = "comm_name") %>%
  as.data.frame()

# Get the sample information
samples <- enframe(names(comms_df)) %>%
  select(value) %>%
  mutate(sMAG=value) %>%
  column_to_rownames(var = "value") %>%
  as.data.frame()

# create a phyloseq object
# otu_table()   OTU Table:         [ 422491 taxa and 713 samples ]
# sample_data() Sample Data:       [ 713 samples by 1 sample variables ]
comms_phyloseq <- phyloseq(otu_table(comms_df, taxa_are_rows = TRUE), sample_data(samples))

# We will filter the tables using differen thresholds of the GC prevalence
prevs <- c(0.005, seq(0.01, 0.05, by = 0.01))
res <- future_map(prevs, function(X){
  filter_prev(physeq = comms_phyloseq, prev = X, data = comms)
}, .progress = TRUE
)

names(res) <- as.character(prevs)

# Let's write it to disk
dev_null <- future_map(res, function(X){
  prev <- X$mag_summary$prev %>% unique()

  fname <- file.path("results", "tables", paste0("euk_communities_prev-", prev, ".tsv.gz"))
  write_tsv(X$data, file = fname, col_names = TRUE)

  fname1 <- file.path("results", "tables", paste0("euk_communities_prev-", prev, "_gene-diff.tsv.gz"))
  write_tsv(X$mag_summary, file = fname1, col_names = TRUE)
}, .progress = TRUE)


# What we want to do?
#  1) get all genes/comms that belong to an eggnog in the known space
#  2) filter all by a certain prevalence
#  3) create tables

comm_k_og <- all_cls_annot %>%
  filter(grepl("K", category)) %>%
  mutate(cl_name = as.character(cl_name))

# We find consensus eggNOG annotations for each gene cluster community
comm_k_og_mv <- get_majority(comm_k_og %>% select(comm_name,bestOG) %>% mutate(bestOG = ifelse(is.na(bestOG), "NONE", bestOG))) %>%
  rename(bestOG_mv = majority)

all_cls_annot_mv <- all_cls_annot %>%
  left_join(comm_k_og_mv %>% mutate(bestOG_mv = ifelse(bestOG_mv == "NONE", NA, bestOG_mv)))

all_cls_annot_mv %>%
  filter(category == "GU", COG_categories == "S")


comms_e_n <- all_cls_annot_mv %>%
  mutate(bestOG_mv = ifelse(is.na(bestOG_mv), comm_name, bestOG_mv),
         bestOG_mv = ifelse((is.na(bestOG_mv) & category == "SINGL"), bestOG, bestOG_mv)) %>%
  filter(category != "DISC", !is.na(bestOG_mv)) %>%
  select(sMAG,  bestOG_mv) %>%
  group_by(sMAG,  bestOG_mv) %>%
  count() %>%
  ungroup()


# Let's explore their prevalence
comms_e_df <- comms_e_n %>%
  spread(key = sMAG, value = n, fill = 0) %>%
  column_to_rownames(var = "bestOG_mv") %>%
  as.data.frame()

samples <- enframe(names(comms_e_df)) %>%
  select(value) %>%
  mutate(sMAG=value) %>%
  column_to_rownames(var = "value") %>%
  as.data.frame()

comms_e_phyloseq <- phyloseq(otu_table(comms_e_df, taxa_are_rows = TRUE), sample_data(samples))


prevs <- c(0.005, seq(0.01, 0.05, by = 0.01))
res <- future_map(prevs, function(X){
  filter_prev_og(physeq = comms_e_phyloseq, prev = X, data = comms)
}, .progress = TRUE
)
names(res) <- as.character(prevs)

# Let's write it to disk
dev_null <- future_map(res, function(X){
  prev <- X$mag_summary$prev %>% unique()

  fname <- file.path("results", "tables", paste0("euk_communities-OG_prev-", prev, ".tsv.gz"))
  write_tsv(X$data, file = fname, col_names = TRUE)

  fname1 <- file.path("results", "tables", paste0("euk_communities-OG_prev-", prev, "_gene-diff.tsv.gz"))
  write_tsv(X$mag_summary, file = fname1, col_names = TRUE)
}, .progress = TRUE)



