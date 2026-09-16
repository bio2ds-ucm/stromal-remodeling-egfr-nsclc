# DESCRIPTION ----
# Osimertinib resistance cell line models.
# Cell lines HCC827 and H1975
# Normalization and rlog transformation of RNA-sequencing data.

# LIBRARIES ----
library(readxl)
library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)

# HCC827 ----

# 1. Load raw count matrix ----

# Read .xlsx
raw_hcc827 <- read_excel("./Data/Cell_line_models/Osimertinib_resistance/HCC827_rawcounts.xlsx")

# Set gene names as rownames
gene <- raw_hcc827$Gene

raw_hcc827 <- raw_hcc827 |>
  select(-Gene)

rownames(raw_hcc827) <- gene

# 2. Sensitivity/resistance groups ----

# Add column for clone id
metadata_hcc827 <- tibble(clone = colnames(raw_hcc827))

# Sensitivity/resistance groups
metadata_hcc827 <- metadata_hcc827 |>
  mutate(
    parental = case_when(
      str_detect(clone, "MC2") ~ "MC1 (monoclonal)",
      str_detect(clone, "MC3") ~ "MC2 (monoclonal)",
      str_detect(clone, "MC6") ~ "MC3 (monoclonal)"
    ),
    sensitivity = case_when(str_detect(clone, "par") ~ "sensitive", .default = "resistant") |>
      factor() |>
      relevel(ref = "sensitive")
  )

# Add clone ids as rownames
rownames(metadata_hcc827) <- metadata_hcc827$clone

# 3. Rlog transformation ----

# Filter genes with raw counts = 0 in all samples
# Number of samples with counts > 0 (per gene)
nz_genes_per_sample <- raw_hcc827 |>
  apply(
    MARGIN = 1,
    FUN = function(x)
      sum(x > 0)
  )

# Genes with counts > 0 in at least one sample
nz_genes <- names(nz_genes_per_sample[nz_genes_per_sample > 0])

# Check dimensions
raw_hcc827 |>
  dim()

# Number of genes with counts > 0 in at least one sample
nz_genes |>
  length()

# Filter out genes with 0 counts in all samples
raw_hcc827 <- raw_hcc827[nz_genes, metadata_hcc827$clone]
rownames(raw_hcc827) <- nz_genes

# DESeq object

# Check dimensions
raw_hcc827 |>
  dim()

metadata_hcc827 |>
  dim()

# Check sample order
raw_hcc827 <- raw_hcc827[, rownames(metadata_hcc827)]

all(colnames(raw_hcc827) == rownames(metadata_hcc827))

# Set gene names as rownames
rownames(raw_hcc827) <- nz_genes

# rlog transformation
dss_hcc827 <- DESeqDataSetFromMatrix(countData = raw_hcc827,
                                     colData = metadata_hcc827,
                                     design = ~ sensitivity)

dss_hcc827 <- estimateSizeFactors(dss_hcc827)
norm_rlog_hcc827 <- rlog(dss_hcc827, blind = FALSE)

# Merge rlog-transformed counts and metadata ----

norm_rlog_hcc827_metadata <- metadata_hcc827 |>
  merge(
    assay(norm_rlog_hcc827) |>
      t() |>
      as.data.frame() |>
      rownames_to_column(var = "clone"),
    by = "clone"
  )

# Clear environment
rm(list = setdiff(ls(), "norm_rlog_hcc827_metadata"))

# H1975 ----

# 1. Load raw count matrix ----

# Read .xlsx
raw_h1975 <- read_excel("./Data/Cell_line_models/Osimertinib_resistance/H1975_rawcounts.xlsx")
raw_h1975 |>
  View()

# Set gene names as rownames
gene <- raw_h1975$Gene

raw_h1975 <- raw_h1975 |>
  select(-Gene)

rownames(raw_h1975) <- gene

# 2. Sensitivity/resistance groups ----

# Add column for clone id
metadata_h1975 <- tibble(clone = colnames(raw_h1975))

# Sensitivity/resistance groups
metadata_h1975 <- metadata_h1975 |>
  mutate(sensitivity = ifelse(str_detect(clone, "wt"), "sensitive", "resistant"))

# Add clone ids as rownames
rownames(metadata_h1975) <- metadata_h1975$clone

# 3. Rlog transformation ----

# Prior to normalization, filter genes with raw counts = 0 in all samples
# Number of samples with counts > 0 (per gene)
nz_genes_per_sample <- raw_h1975 |>
  apply(
    MARGIN = 1,
    FUN = function(x)
      sum(x > 0)
  )

# Genes with counts > 0 in at least one sample
nz_genes <- names(nz_genes_per_sample[nz_genes_per_sample > 0])

# Check dimensions
raw_h1975 |>
  dim()

# Number of genes with counts > 0 in at least one sample
nz_genes |>
  length()

# Filter out genes with 0 counts in all samples
raw_h1975 <- raw_h1975[nz_genes, metadata_h1975$clone]
rownames(raw_h1975) <- nz_genes

# DESeq object

# Check dimensions
raw_h1975 |>
  dim()

metadata_h1975 |>
  dim()

# Check sample order
raw_h1975 <- raw_h1975[, rownames(metadata_h1975)]
all(colnames(raw_h1975) == rownames(metadata_h1975))

# Gene names as rownames
rownames(raw_h1975) <- nz_genes

# rlog transformation
dss_h1975 <- DESeqDataSetFromMatrix(countData = raw_h1975,
                                    colData = metadata_h1975,
                                    design = ~ sensitivity)

dss_h1975 <- estimateSizeFactors(dss_h1975)
norm_rlog_h1975 <- rlog(dss_h1975, blind = FALSE)

# Merge rlog-transformed counts and metadata ----

norm_rlog_h1975_metadata <- metadata_h1975 |>
  merge(
    assay(norm_rlog_h1975) |>
      t() |>
      as.data.frame() |>
      rownames_to_column(var = "clone"),
    by = "clone"
  )

# Clear environment
rm(list = setdiff(
  ls(),
  c("norm_rlog_hcc827_metadata", "norm_rlog_h1975_metadata")
))

# SAVE ----

# HCC827
saveRDS(
  norm_rlog_hcc827_metadata,
  "./Data/Cell_line_models/Osimertinib_resistance/HCC827_normrlog_counts.rds"
)

# H1975
saveRDS(
  norm_rlog_h1975_metadata,
  "./Data/Cell_line_models/Osimertinib_resistance/H1975_normrlog_counts.rds"
)
