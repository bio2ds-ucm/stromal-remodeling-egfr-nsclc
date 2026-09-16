# DESCRIPTION ----
# OSIRESP cohort, including pre-treatment biopsies AND additionally collected
# post-progression biopsies.
# Analyses of post-progression versus pre-treatment paired samples by compartment
# (i.e., tumor or stroma).

# LIBRARIES ----
library(openxlsx)
library(GeomxTools)
library(dplyr)
library(msigdbr)
library(limma)
library(clusterProfiler)

# CTA GENE SETS ----

# Load Bruker CTA gene sets
cta_gs <- readRDS("./Data/Bruker_CTA_assay_gene_sets.rds")

# Remove gene sets corresponding to organ of origin
cta_gs <- cta_gs |>
  filter(!(
    gs_name %in% c("Glioma", "Leukemia", "Melanoma", "Prostate Cancer")
  ))

# YAP CONSERVEED SIGNATURE ----

# Documentation:
# https://www.gsea-msigdb.org/gsea/msigdb/cards/CORDENONSI_YAP_CONSERVED_SIGNATURE.html

msig <- msigdbr(species = "Homo sapiens")

# Get genes
yap_genes <- msig |>
  filter(gs_name == "CORDENONSI_YAP_CONSERVED_SIGNATURE") |>
  pull(gene_symbol) |>
  unique()

# Check set size
yap_genes |>
  length()

# Add to CTA gene sets
yap_taz_df <- tibble(gene_symbol = yap_genes) |>
  mutate(gs_name = "Cordenonsi YAP Conserved Signature")

cta_gs <- list(cta_gs, yap_taz_df) |>
  bind_rows()

# DATA LOADING ----

# Load NanoStringGeoMxSet object
osiresp_rb <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_q3norm.rds")

# Extract metadata from NanoStringGeoMxSet object
osiresp_rb <- pData(osiresp_rb)

# Expression data + biopsy type annotation (pre-treatment or at progression)
expr_data <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_vst_medianvalues.rds")

# DATA PREPARATION ----

# Order levels of biopsy_stye
expr_data <- expr_data |>
  mutate(biopsy_type = biopsy_type |>
           factor(
             levels = c("pre-treatment", "at progression"),
             labels = c("pre", "post")
           )) |>
  mutate(patient_biopsy = paste(patient_id, biopsy_type, sep = "_"))

# Expression data by compartment
expr_stroma <- expr_data |>
  filter(compartment == "Stroma")

expr_tumor <- expr_data |>
  filter(compartment == "Tumor")

# Patients with paired samples (stroma compartment)
n_biopsies_stroma <- expr_stroma |>
  group_by(patient_id) |>
  summarize(n_biopsies = n()) |>
  filter(n_biopsies == 2)

# Patients with paired samples (tumor compartment)
n_biopsies_tumor <- expr_tumor  |>
  group_by(patient_id) |>
  summarize(n_biopsies = n()) |>
  filter(n_biopsies == 2)

# DIFFERENTIAL EXPRESSION ANALYSIS ----

# Limma-trend workflow.
# Only paired samples per compartment considered.

# 1. Stroma compartment ----

# Subset paired data
paired_stroma_pat <- n_biopsies_stroma$patient_id

expr_paired_stroma <- expr_stroma |>
  filter(patient_id %in% paired_stroma_pat) |>
  arrange(patient_id, biopsy_type)

osiresp_rb_paired_stroma <- osiresp_rb |>
  filter(patient_id %in% paired_stroma_pat) |>
  select(patient_id) |>
  unique() |>
  inner_join(expr_paired_stroma |>
               select(patient_id, biopsy_type)) |>
  arrange(patient_id, biopsy_type)

# Check order
all(expr_paired_stroma$patient_id == osiresp_rb_paired_stroma$patient_id)
all(expr_paired_stroma$biopsy_type == osiresp_rb_paired_stroma$biopsy_type)

# Design matrix
design_stroma <- model.matrix( ~ patient_id + biopsy_type, data = osiresp_rb_paired_stroma)

# Expression matrix
expr_mat_stroma <- expr_paired_stroma |>
  ungroup() |>
  select(-patient_id, -compartment, -biopsy_type, -patient_biopsy) |>
  t()

# Set colnames
colnames(expr_mat_stroma) <- expr_paired_stroma$patient_biopsy

# Fit
fit_stroma <- lmFit(expr_mat_stroma, design_stroma, robust = TRUE)
fit_stroma <- eBayes(fit_stroma, trend = TRUE, robust = TRUE)

# Results
fit_stroma_res <- topTable(fit_stroma,
                           coef = "biopsy_typepost",
                           number = Inf,
                           sort.by = "P")

# 2. Tumor compartment ----

# Subset paired data
paired_tumor_pat <- n_biopsies_tumor$patient_id

expr_paired_tumor <- expr_tumor |>
  filter(patient_id %in% paired_tumor_pat) |>
  arrange(patient_id, biopsy_type)

osiresp_rb_paired_tumor <- osiresp_rb |>
  filter(patient_id %in% paired_tumor_pat) |>
  select(patient_id) |>
  unique() |>
  inner_join(expr_paired_tumor |>
               select(patient_id, biopsy_type)) |>
  arrange(patient_id, biopsy_type)

# Check order
all(expr_paired_tumor$patient_id == osiresp_rb_paired_tumor$patient_id)
all(expr_paired_tumor$biopsy_type == osiresp_rb_paired_tumor$biopsy_type)

# Design matrix
design_tumor <- model.matrix( ~ patient_id + biopsy_type, data = osiresp_rb_paired_tumor)

# Expression matrix
expr_mat_tumor <- expr_paired_tumor |>
  ungroup() |>
  select(-patient_id, -compartment, -biopsy_type, -patient_biopsy) |>
  t()

# Set colnames
colnames(expr_mat_tumor) <- expr_paired_tumor$patient_biopsy

# Fit
fit_tumor <- lmFit(expr_mat_tumor, design_tumor, robust = TRUE)
fit_tumor <- eBayes(fit_tumor, trend = TRUE, robust = TRUE)

# Results
fit_tumor_res <- topTable(fit_tumor,
                          coef = "biopsy_typepost",
                          number = Inf,
                          sort.by = "P")

# GSEA ----

# Seed for replicability
seed <- 1850

# 1. Stroma compartment ----

# GSEA ranking metric
dgea_metric_stroma <- fit_stroma_res$t
names(dgea_metric_stroma) <- rownames(fit_stroma_res)
dgea_metric_stroma <- sort(dgea_metric_stroma, decreasing = TRUE)

# Set seed for replicability
set.seed(seed)
gsea_stroma <- GSEA(
  geneList = dgea_metric_stroma,
  TERM2GENE = cta_gs,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

# 2. Tumor compartment ----

# GSEA ranking metric
dgea_metric_tumor <- fit_tumor_res$t
names(dgea_metric_tumor) <- rownames(fit_tumor_res)
dgea_metric_tumor <- sort(dgea_metric_tumor, decreasing = TRUE)

# Set seed for replicability
set.seed(seed)
gsea_tumor <- GSEA(
  geneList = dgea_metric_tumor,
  TERM2GENE = cta_gs,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

# SUPP FILE WITH RESULTS ----

dgea_gsea_res <- list(
  `Stroma (paired DGEA limma)` = fit_stroma_res |> arrange(desc(logFC)),
  `Tumor (paired DGEA limma)` = fit_tumor_res |> arrange(desc(logFC)),
  `Stroma (paired GSEA)` = gsea_stroma@result,
  `Tumor (paired = GSEA)` = gsea_tumor@result
)

write.xlsx(dgea_gsea_res, rowNames = TRUE, file = "./Results/Tables/9_OSIRESP_RB_postprogression_vs_pretreatment.xlsx")

# SAVE RESULTS ----

# Clear environment
rm(list = setdiff(
  ls(),
  c(
    "expr_paired_stroma",
    "expr_paired_tumor",
    "fit_stroma_res",
    "fit_tumor_res",
    "gsea_stroma",
    "gsea_tumor"
  )
))

# save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/9_OSIRESP_RB_postprogression_vs_pretreatment.RData")
save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/9_OSIRESP_RB_postprogression_vs_pretreatment.RData")
