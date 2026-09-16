# DESCRIPTION ----
# Osimertinib drug-tolerant persistence cell line models from
# Criscione et al. (https://pubmed.ncbi.nlm.nih.gov/36575215/)
# RNA-sequencing data publicly available at NCBI's GEO (accession: GSE193258)

# HCC827 line: initial processing of RNA-sequencing data and metadata
library(GEOquery)
library(dplyr)
library(stringr)

# LOAD DATA FROM GEO ----

geo_data <- getGEO("GSE193258", GSEMatrix = TRUE)

geo_data |>
  View()

# 1. Metadata ----

# Metadata
metadata <- phenoData(geo_data[[1]])@data

# Explore metadata
metadata |>
  class()

metadata |>
  dim()

metadata |>
  View()

metadata$`cell line:ch1` |>
  table(useNA = "always")

# Subset metadata for HCC827
metadata <- metadata |>
  filter(`cell line:ch1` == "HCC827")

metadata <- metadata |>
  filter(str_detect(title, "DMSO|DTP"))

metadata <- metadata |>
  dplyr::select(title, `treatment:ch1`) |>
  dplyr::rename(treatment = `treatment:ch1`)

# 2. Expression data ----

# Raw count matrix (.tsv file downloaded from GEO)
expr_raw <- read.table(file = "./Data/Cell_line_models/Persistence/GSE193258_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE)

expr_raw |>
  dim()

# Check sample IDs correspondence between expression data and metadata
all(rownames(metadata) %in% colnames(expr_raw))
setdiff(rownames(metadata), colnames(expr_raw))
metadata[setdiff(rownames(metadata), colnames(expr_raw)), ]

# Save gene IDs in vector
gene_id <- expr_raw$GeneID

# Subset expression data (only samples of interest)
expr_raw <- expr_raw[, intersect(rownames(metadata), colnames(expr_raw))]

# Set gene IDs as rownames
rownames(expr_raw) <- gene_id

# Check dimensions
metadata |>
  dim()

expr_raw |>
  dim()

# FILTER GENES WITH ZERO COUNTS IN ALL SAMPLES ----

# Filter genes with zero counts in all samples
# Number of samples with counts > 0 (per gene)
non_zero_n <- expr_raw |>
  apply(
    MARGIN = 1,
    FUN = function(x)
      sum(x > 0)
  )

# Genes with > 0 counts in at least one sample
non_zero_subset <- names(non_zero_n[non_zero_n > 0])

non_zero_subset |>
  length()

expr_raw |>
  dim()

# Filter out genes with 0 counts in all samples
expr_raw <- expr_raw[non_zero_subset, ]

expr_raw |>
  dim()

# SAVE METADATA AND EXPR DATA ----

list_data <- list(metadata = metadata, expr_raw = expr_raw)

saveRDS(list_data,
        "./Data/Cell_line_models/Persistence/GSE193258_HCC827.rds")
