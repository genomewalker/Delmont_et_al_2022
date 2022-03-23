library(tidyverse)
library(RSQLite)
library(maditr)
library(tidytext)
library(furrr)
library(doFuture)
options(future.globals.maxSize=10485760000)

registerDoFuture()
# You might want to change it
ncores <- 4
plan(future::multisession, workers = ncores)

source("libs/lib.R")

# Get it from http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
# Change the version to the version you used
eggdb <- "data/eggnog-db/eggnog.db"
con <- RSQLite::dbConnect(RSQLite::SQLite(), eggdb)

# Get the results from the eggnog-mapper
results_dir <- "data/eggnog-annotations/01_Results_EggNog_annotationV1"
emapper_files_pattern <-
  ".metaeuk.eggnog.annotation.emapper.annotations"
emapper_files_pattern_filt <-
  paste0(emapper_files_pattern, "-filt.tsv")
emapper_files <-
  list.files(path = results_dir,
             pattern = emapper_files_pattern,
             full.names = TRUE)

# Get taxonomic information
taxmap <- read_tsv(
  "data/eggnog-db/emapper-v2.taxmap.tsv",
  col_names = c("taxid", "tax_name"),
  col_types = list("c", "c"),
  comment = "#"
)

# Get the COG funcats
funcat <- fread("data/eggnog-db/eggnog-funcat.tsv", ) %>%
  setNames(c("COG_Functional_Category", "cog_category", "cog_description"))


# Get annotations
annotations <- tbl(con, "og") %>% collect()



# Read and process the emmaper results - RAW
emapper_results <- future_map_dfr(emapper_files,
             read_emmaper_files,
             pattern = emapper_files_pattern,
             .progress = TRUE)

# Read and process the emmaper results - LOWEST OG
emapper_results_lower <- future_map_dfr(emapper_files,
             read_emmaper_files_lower,
             pattern = emapper_files_pattern,
             .progress = TRUE)

# Get SMAGs
smags <- emapper_results$sMAG %>% unique()
smags_lower <- emapper_results_lower$sMAG %>% unique()


# Save results
dev_null <- future_map(smags_lower, function(X) {
  fname <-
    file.path(
      "results", "lowestOG",
      paste0(X, emapper_files_pattern_filt, ".gz")
    )
  data <- emapper_results_lower %>%
    write_tsv(data, file = fname)
}, .progress = TRUE)

# Let's get the names of each eggnog
emapper_results_descs <- emapper_results %>%
  select(eggNOG_free_text_description, description) %>%
  distinct()

emapper_results_descs %>%
  filter(eggNOG_free_text_description != description) %>%
  nrow()

text <- emapper_results_lower %>%
  select(description) %>%
  distinct() %>%
  mutate(text = gsub("_", " ", description)) %>%
  unnest_tokens(lines, text, token = "lines", drop = TRUE)

# Aggregate results and save to DB/file
emapper_results_agg <- emapper_results_lower %>%
  inner_join(text) %>%
  rename(description_clean = lines) %>%
  select(sMAG, description_clean) %>%
  group_by(sMAG, description_clean) %>%
  count(sort = T) %>%
  ungroup()

write_tsv(
  emapper_results_agg,
  "results/lowestOG/01_Results_EggNog_annotationV1-filt_agg-clean.tsv.gz"
)



