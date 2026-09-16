# DESCRIPTION ----
# Analyses of osimertinib drug-tolerant persistence cell line models from
# Criscione et al. (https://pubmed.ncbi.nlm.nih.gov/36575215/)
# RNA-sequencing data publicly available at NCBI's GEO (accession: GSE193258)

# H1975, HCC827, and HCC2935 cell lines (independently analyzed)

# - Differential Gene Expression Analysis (DGEA), limma-voom approach
# - Gene Set Enrichment Analysis (GSEA) using Reactome pathways

# LIBRARIES ----
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(edgeR)
library(clusterProfiler)

# REACTOME GENE SETS ----

# Reactome gene sets
gs_react <- msigdbr(species = "Homo sapiens",
                    collection  = "C2",
                    subcollection = "CP:REACTOME")

msig <- msigdbr(species = "Homo sapiens")

# Add Cordenonsi YAP Conserved Signature
yap_taz_genes <- msig |>
  filter(gs_name == "CORDENONSI_YAP_CONSERVED_SIGNATURE") |>
  pull(gene_symbol) |>
  unique()

yap_taz_df <- tibble(gene_symbol = yap_taz_genes) |>
  mutate(gs_name = "CORDENONSI_YAP_CONSERVED_SIGNATURE")

gs_react <- list(gs_react, yap_taz_df) |>
  bind_rows()

# Gene annotation
conv <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gs_react$gene_symbol),
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

gs_react <- merge(gs_react, conv, by.x = "gene_symbol", by.y = "SYMBOL")

# Select columns
gs_react <- gs_react |>
  dplyr::select(gs_name, ENTREZID)

# Filter out gene sets with size < 15
gs_react_count <- gs_react |>
  group_by(gs_name) |>
  summarise(count = n())

min_count <- 15
gs_react_count <- gs_react_count |>
  filter(count >= min_count)

gs_react <- gs_react |>
  filter(gs_name %in% gs_react_count$gs_name)

# Clear environment
rm(list = setdiff(ls(), "gs_react"))

# H1975 cell line ----

# Load metadata and expression data
data <- readRDS("./Data/Cell_line_models/Persistence/GSE193258_H1975.rds")
metadata <- data$metadata
expr <- data$expr_raw

# 1. DGEA (limma-voom) ----

# Design matrix
design <- model.matrix(~ 1 + treatment, data = metadata)

design |>
  head()

# Create DGEList and normalize
y <- DGEList(expr)
y <- calcNormFactors(y)

# Voom transformation
v <- voom(y, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design, robust = TRUE)
fit <- eBayes(fit, robust = TRUE)

fit_results <- topTable(fit, n = Inf, p = 1)

# 2. GSEA ----

# Ranking metric
dgea_metric <- fit_results$t
names(dgea_metric) <- rownames(fit_results)
dgea_metric <- sort(dgea_metric, decreasing = TRUE)

# GSEA
set.seed(1850)
gsea <- GSEA(
  geneList = dgea_metric,
  TERM2GENE = gs_react,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0,
  nPermSimple = 100000
)

gsea_h1975 <- gsea
fit_results_h1975 <- fit_results

# Clear environment
rm(list = setdiff(ls(), c(
  "gs_react", "gsea_h1975", "fit_results_h1975"
)))

# HCC827 cell line ----

# Load metadata and expression data
data <- readRDS("./Data/Cell_line_models/Persistence/GSE193258_HCC827.rds")

metadata <- data$metadata
expr <- data$expr_raw

# 1. DGEA (limma-voom) ----

# Design matrix
design <- model.matrix(~ 1 + treatment, data = metadata)

design |>
  head()

# Create DGEList and normalize
y <- DGEList(expr)
y <- calcNormFactors(y)

# Voom transformation
v <- voom(y, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design, robust = TRUE)
fit <- eBayes(fit, robust = TRUE)

fit_results <- topTable(fit, n = Inf, p = 1)

# 2. GSEA ----

# Ranking metric
dgea_metric <- fit_results$t
names(dgea_metric) <- rownames(fit_results)
dgea_metric <- sort(dgea_metric, decreasing = TRUE)

# GSEA
set.seed(1850)
gsea <- GSEA(
  geneList = dgea_metric,
  TERM2GENE = gs_react,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0,
  nPermSimple = 100000
)

gsea_hcc827 <- gsea
fit_results_hcc827 <- fit_results

# Clear environment
rm(list = setdiff(
  ls(),
  c(
    "gs_react",
    "gsea_h1975",
    "fit_results_h1975",
    "gsea_hcc827",
    "fit_results_hcc827"
  )
))

# HCC2935 cell line ----

# Load metadata and expression data
data <- readRDS("./Data/Cell_line_models/Persistence/GSE193258_HCC2935.rds")

metadata <- data$metadata
expr <- data$expr_raw

# 1. DGEA (limma-voom) ----

# Design matrix
design <- model.matrix(~ 1 + treatment, data = metadata)

design |>
  head()

# Create DGEList and normalize
y <- DGEList(expr)
y <- calcNormFactors(y)

# Voom transformation
v <- voom(y, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design, robust = TRUE)
fit <- eBayes(fit, robust = TRUE)

fit_results <- topTable(fit, n = Inf, p = 1)

# 2. GSEA ----

# Ranking metric
dgea_metric <- fit_results$t
names(dgea_metric) <- rownames(fit_results)
dgea_metric <- sort(dgea_metric, decreasing = TRUE)

# GSEA
set.seed(1850)
gsea <- GSEA(
  geneList = dgea_metric,
  TERM2GENE = gs_react,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0,
  nPermSimple = 100000
)

gsea_hcc2935 <- gsea
fit_results_hcc2935 <- fit_results

# Clear environment
rm(list = setdiff(
  ls(),
  c(
    "gsea_h1975",
    "fit_results_h1975",
    "gsea_hcc827",
    "fit_results_hcc827",
    "gsea_hcc2935",
    "fit_results_hcc2935"
  )
))

# SAVE RESULTS ----

save.image(
  "./Results/Intermediate/Cell_line_models/Persistence/1_GSE193258_cell_lines_GSEA.RData"
)
