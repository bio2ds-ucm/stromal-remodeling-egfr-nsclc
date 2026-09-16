# DESCRIPTION ----
# Analysis of osimertinib resistance cell line models.
# Cell lines HCC827 and H1975 (independently analyzed)

# HCC827 cell line:
# - Three parental monoclonal cell lines, each derived three resistant clones
# Analyses for each parental monoclonal cell line performed independently

# H1975:
# - Parental polyclonal cell line deriving resistant clones

# Gene Set Enrichment Analysis (GSEA) using Reactome pathways

# LIBRARIES ----
library(dplyr)
library(tidyr)
library(msigdbr)
library(clusterProfiler)

# REACTOME GENE SETS ----

# Reactome gene sets
gs_react <- msigdbr(species = "Homo sapiens",
                    collection  = "C2",
                    subcollection = "CP:REACTOME")

# Select columns
gs_react <- gs_react |>
  select(gs_name, gene_symbol)

# Add Cordenonsi YAP Conserved Signature
msig <- msigdbr(species = "Homo sapiens")

yap_taz_genes <- msig |>
  filter(gs_name == "CORDENONSI_YAP_CONSERVED_SIGNATURE") |>
  pull(gene_symbol) |>
  unique()

yap_taz_df <- tibble(gene_symbol = yap_taz_genes) |>
  mutate(gs_name = "CORDENONSI_YAP_CONSERVED_SIGNATURE")

gs_react <- list(gs_react, yap_taz_df) |>
  bind_rows()

# Number of genes per gene set
gs_react_count <- gs_react |>
  group_by(gs_name) |>
  summarise(count = n())

# Filter out gene sets with size < 15
min_count <- 15
gs_react_count <- gs_react_count |>
  filter(count >= min_count)

gs_react <- gs_react |>
  filter(gs_name %in% gs_react_count$gs_name)

# Clear environment
rm(list = setdiff(ls(), c("gs_react")))

# HCC827 ----

# 1. Load rlog-transformed count matrix ----

expr_hcc827 <- readRDS("./Data/Cell_line_models/Osimertinib_resistance/HCC827_normrlog_counts.rds")

# 2. Compute log2FC ----

# log2FC computed as:
# (mean expression in resistant clones) - (expression in parental cell line)
# after rlog transformation

# Data in long format
expr_hcc827_long <- expr_hcc827 |>
  pivot_longer(
    cols = -c(clone, parental, sensitivity),
    names_to = "gene",
    values_to = "value"
  )

# Mean value for resistant clones (independenttly for each parental cell line)
clone_hcc827 <- expr_hcc827_long |>
  filter(sensitivity == "resistant") |>
  group_by(parental, gene) |>
  summarise(mean_clone = mean(value))

# Get expression values in parental cell lines
parental_hcc827 <- expr_hcc827_long |>
  filter(sensitivity == "sensitive") |>
  select(parental, gene, parental_value = value)

# Join and compute log2FC
log2fc_hcc827 <- parental_hcc827 |>
  left_join(clone_hcc827, by = c("parental", "gene")) |>
  mutate(log2fc = mean_clone - parental_value)

# Split dataset to get one per parental cell line
log2fc_hcc827 <- log2fc_hcc827 |>
  group_split(parental)

# Set names for each parental monoclonal cell line
names(log2fc_hcc827) <- c("MC1", "MC2", "MC3")

# Clear environment
rm(list = setdiff(ls(), c("gs_react", "log2fc_hcc827")))

# 3. GSEA ----

# Seed for replicability
seed <- 1850

# Create empty list (to save GSEA results)
gsea_hcc827 <- vector(mode = "list", length = length(log2fc_hcc827))
names(gsea_hcc827) <- names(log2fc_hcc827)

for (ind in 1:length(log2fc_hcc827)) {
  # Ranking metric
  dgea_metric_hcc827 <- log2fc_hcc827[[ind]]$log2fc
  names(dgea_metric_hcc827) <- log2fc_hcc827[[ind]]$gene
  dgea_metric_hcc827 <- sort(dgea_metric_hcc827, decreasing = TRUE)
  
  # GSEA: Reactome gene sets
  # Set seed for replicability
  set.seed(seed)
  gsea <- GSEA(
    geneList = dgea_metric_hcc827,
    TERM2GENE = gs_react,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0,
    nPermSimple = 100000
  )
  
  gsea_hcc827[[ind]] <- gsea
  
}

# Clear environment
rm(list = setdiff(ls(), c(
  "gs_react", "log2fc_hcc827", "gsea_hcc827"
)))

# H1975 ----

# 1. Load rlog-transformed count matrix ----

rlog_h1975 <- readRDS("./Data/Cell_line_models/Osimertinib_resistance/H1975_normrlog_counts.rds")

# 2. Compute log2FC ----

# log2FC computed as:
# (mean expression in resistant clones) - (expression in parental cell line)
# after rlog transformation

# Data in long format
rlog_h1975_long <- rlog_h1975 |>
  pivot_longer(
    cols = -c(clone, sensitivity),
    names_to = "gene",
    values_to = "value"
  )

# Mean value for resistant clones
clone_h1975 <- rlog_h1975_long |>
  filter(sensitivity == "resistant") |>
  group_by(gene) |>
  summarise(mean_clone = mean(value))

# Get expression values in parental cell line
parental_h1975 <- rlog_h1975_long |>
  filter(sensitivity == "sensitive") |>
  select(gene, parental_value = value)

# Join and compute log2fc
log2fc_h1975 <- parental_h1975 |>
  left_join(clone_h1975, by = "gene") |>
  mutate(log2fc = mean_clone - parental_value)

# Clear environment
rm(list = setdiff(
  ls(),
  c("gs_react", "log2fc_hcc827", "gsea_hcc827", "log2fc_h1975")
))
# 3. GSEA ----

# Seed for replicability
seed <- 1850

# Ranking metric
dgea_metric_h1975 <- log2fc_h1975$log2fc
names(dgea_metric_h1975) <- log2fc_h1975$gene
dgea_metric_h1975 <- sort(dgea_metric_h1975, decreasing = TRUE)

# GSEA: Reactome gene sets
# Set seed for replicability
set.seed(seed)
gsea_h1975 <- GSEA(
  geneList = dgea_metric_h1975,
  TERM2GENE = gs_react,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0,
  nPermSimple = 100000
)

# Clear environment
rm(list = setdiff(
  ls(),
  c("log2fc_hcc827", "gsea_hcc827", "log2fc_h1975", "gsea_h1975")
))

# SAVE ----
save.image(file = "./Results/Intermediate/Cell_line_models/Osimertinib_resistance/1_Osimertinib_resistance_cell_lines_GSEA.RData")
