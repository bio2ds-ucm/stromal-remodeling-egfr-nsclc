# DESCRIPTION -----
# OSIRESP cohort.
# Unsupervised analysis (Leiden clustering) for identification of molecular subtypes
# based on stromal features.

# LIBRARIES ----
library(GeomxTools)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Seurat)
library(uwot)
library(janitor)
library(limma)
library(clusterProfiler)
library(GSVA)
library(openxlsx)
library(msigdbr)

# CTA GENE SETS ----

# Load Bruker CTA gene sets
cta_gs <- readRDS("./Data/Bruker_CTA_assay_gene_sets.rds")

# Supp File with Bruker CTA gene sets
write.xlsx(cta_gs, rowNames = FALSE, file = "./Results/Tables/Bruker_CTA_gene_sets.xlsx")

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
expr_data <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_vst_medianvalues.rds")

pheno_data <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_q3norm_pData.rds")

# DATA PREPARATION ----

# Subset expression values from stroma compartment
expr_stroma <- expr_data |>
  ungroup() |>
  filter(compartment == "Stroma") |>
  select(-compartment, -`Negative Probe`)

# Extract patient ids
stroma_pat_id <- expr_stroma |>
  select(patient_id) |>
  pull()

# Remove patient id column and transpose
expr_stroma <- expr_stroma |>
  select(-patient_id) |>
  t()

# Check dimensions
expr_stroma |>
  dim()

# Set patient ids as colnames
colnames(expr_stroma) <- stroma_pat_id

# MOST VARIABLE GENES ----

# Create seurat object
seurat_obj <- CreateSeuratObject(counts = expr_stroma)

# Counts are already normalized and transformed
LayerData(seurat_obj, layer = "data") <- LayerData(seurat_obj, layer = "counts")

# Find top 500 variable genes
nfeatures <- 500
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = nfeatures)

var_features <- VariableFeatures(seurat_obj)

# PCA ----

# Scale data
seurat_obj <- ScaleData(seurat_obj)

# PCA
npc <- 20
seurat_obj <- RunPCA(seurat_obj, features = var_features, npcs = npc)

# UMAP ----
seed <- 1850
set.seed(seed)
seurat_obj  <- RunUMAP(
  seurat_obj,
  dims = 1:npc,
  n.neighbors = 10,
  min.dist = 0.1
)

# LEIDEN CLUSTERING ----

# Parameters
k <- 10
res <- 1.25
seed <- 1850

set.seed(seed)
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "pca",
  dims = 1:npc,
  k.param = k
)

# algorithm = 4 for Leiden clustering
seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = 4)

# UMAP coordinates
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# Check order
all(names(seurat_obj$seurat_clusters) == rownames(umap_coords))

# Leiden clustering results
leid_res <- umap_coords |>
  as.data.frame() |>
  rownames_to_column(var = "patient_id") |>
  inner_join(
    pheno_data |>
      mutate(stage = case_when(
        str_detect(stage_at_baseline, "IV") ~ "IV",
        str_detect(stage_at_baseline, "III") ~ "III",
        .default = "II"
      )) |>
      select(
        patient_id,
        long_response,
        disease_progression_to_osimertinib,
        pfs_time_months,
        stage,
        egfr_mutation_type,
        tobacco_history
      ) |>
      unique()
  ) |>
  mutate(cluster = paste0("C", seurat_obj$seurat_clusters))

# Cluster size
leid_res$cluster |>
  table()

# ASSOCIATION WITH CLINICAL FEATURES ----

# EGFR mutation type ----

# Two-way frequency table
leid_res |>
  tabyl(egfr_mutation_type, cluster) |>
  adorn_totals(where = c("row", "col")) |>
  adorn_percentages(c("row")) |>
  adorn_pct_formatting(digits = 1) |>
  adorn_ns()

# Fisher's test
leid_egfr_mut <- table(leid_res$cluster,
                       leid_res$egfr_mutation_type,
                       dnn = c("cluster", "egfr mutation"))
leid_egfr_mut

fisher.test(leid_egfr_mut)

# Tobacco history ----

# Two-way frequency table
leid_res |>
  tabyl(tobacco_history, cluster) |>
  adorn_totals(where = c("row", "col")) |>
  adorn_percentages(c("row")) |>
  adorn_pct_formatting(digits = 1) |>
  adorn_ns()

# Fisher's test
leid_tobacco <- table(leid_res$cluster,
                      leid_res$tobacco_history,
                      dnn = c("cluster", "tobacco_history"))
leid_tobacco

fisher.test(leid_tobacco)

# TP53 alterations ----

# Point mutation annotations (EGFR and TP53)
tp53_annot <- read.table(
  "./Data/OSIRESP_cohort/Clinical_annotations/OSIRESP_realvariants_EGFR_TP53.csv",
  header = TRUE
)

# Check patients with point mutation annotations
tp53_annot$patient_id %in% leid_res$patient_id |> table()
setdiff(tp53_annot$patient_id, leid_res$patient_id)

leid_res$patient_id %in% tp53_annot$patient_id |> table()

# EGFR mutation detected
tp53_annot$egfr_mut |> table()
# TP53 mutation detected
tp53_annot$tp53_mut |> table()

# Two-way frequency table
table(
  tp53_annot$egfr_mut,
  tp53_annot$tp53_mut,
  dnn = c("EGFR", "TP53"),
  useNA = "always"
)

# Set TP53 mutation annotation to non-available when EGFR mutation is not detected
# (quality criterion)
tp53_annot <- tp53_annot |>
  mutate(tp53_mut = case_when(tp53_mut == "No" &
                                egfr_mut == "No" ~ "Non-available", .default = tp53_mut))

tp53_annot$egfr_mut |> table()
tp53_annot$tp53_mut |> table()

table(
  tp53_annot$egfr_mut,
  tp53_annot$tp53_mut,
  dnn = c("EGFR", "TP53"),
  useNA = "always"
)

# Add TP53 mutation annotation to Leiden clustering results
leid_res <- leid_res |>
  left_join(tp53_annot) |>
  mutate(tp53_mut = case_when(is.na(tp53_mut) ~ "Non-available", .default = tp53_mut))

leid_res$tp53_mut |> table(useNA = "always")

# Two-way frequency table
leid_res |>
  filter(tp53_mut != "Non-available") |>
  tabyl(tp53_mut, cluster) |>
  adorn_totals(where = c("row", "col")) |>
  adorn_percentages(c("row")) |>
  adorn_pct_formatting(digits = 1) |>
  adorn_ns()

# Fisher's test
tp53_eval_pat <- leid_res$tp53_mut != "Non-available"

leid_tp53 <- table(leid_res$cluster[tp53_eval_pat],
                   leid_res$tp53_mut[tp53_eval_pat],
                   dnn = c("cluster", "TP53 mutation"))

fisher.test(leid_tp53)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS ----

# Limma-trend workflow

# 1. Stroma compartment ----

# Design matrix
design_stroma <- model.matrix( ~ 0 + cluster, data = leid_res)

# Fit
fit_stroma <- lmFit(expr_stroma, design_stroma, robust = TRUE)
fit_stroma <- eBayes(fit_stroma, trend = TRUE, robust = TRUE)

# Results
fit_stroma_res <- topTable(fit_stroma, n = Inf, p = 1)

# Supp File with DGEA results
write.xlsx(fit_stroma_res, rowNames = TRUE, file = "./Results/Tables/OSIRESP_stromal_subtypes_dgea_stroma.xlsx")

# 2. Tumor compartment ----

# Subset expression values from tumor compartment
expr_tumor <- expr_data |>
  ungroup() |>
  filter(compartment == "Tumor") |>
  select(-compartment, -`Negative Probe`)

# Check dimensions
expr_tumor |>
  dim()

# Filter patients with available stromal segments
expr_tumor <- expr_tumor |>
  filter(patient_id %in% colnames(expr_stroma))

# Check dimensions
expr_tumor |>
  dim()

# Extract patient ids
tumor_pat_id <- expr_tumor |>
  select(patient_id) |>
  pull()

# Remove patient id column and transpose
expr_tumor <- expr_tumor |>
  select(-patient_id) |>
  t()

# Check dimensions
expr_tumor |>
  dim()

# Set patient ids as colnames
colnames(expr_tumor) <- tumor_pat_id

# Leiden results for patients with expression data from tumor compartment
leid_res_tumor <- leid_res |>
  filter(patient_id %in% tumor_pat_id)

# Check dimensions
leid_res_tumor |>
  dim()

# Design matrix
design_tumor <- model.matrix( ~ 0 + cluster, data = leid_res_tumor)

# Fit
fit_tumor <- lmFit(expr_tumor, design_tumor, robust = TRUE)
fit_tumor <- eBayes(fit_tumor, trend = TRUE, robust = TRUE)

# Results
fit_tumor_res <- topTable(fit_tumor, n = Inf, p = 1)

# Supp File with DGEA results
write.xlsx(fit_tumor_res, rowNames = TRUE, file = "./Results/Tables/OSIRESP_stromal_subtypes_dgea_tumor.xlsx")

# MARKER GENES ----

# Limma-trend workflow
# Contrasts: each cluster vs the others

# 1. Stroma compartment ----

# Design matrix
design_stroma <- model.matrix( ~ 0 + cluster, data = leid_res)

fit_stroma <- lmFit(expr_stroma, design_stroma, robust = TRUE)

# Contrasts: each cluster vs the others
contrasts_mrkr_stroma <- makeContrasts(
  c1_vs_other = clusterC1 - (clusterC2 + clusterC3 + clusterC4) / 3,
  c2_vs_other = clusterC2 - (clusterC1 + clusterC3 + clusterC4) / 3,
  c3_vs_other = clusterC3 - (clusterC1 + clusterC2 + clusterC4) / 3,
  c4_vs_other = clusterC4 - (clusterC1 + clusterC2 + clusterC3) / 3,
  levels = design_stroma
)

fit_stroma <- contrasts.fit(fit_stroma, contrasts_mrkr_stroma)
fit_stroma <- eBayes(fit_stroma, trend = TRUE, robust = TRUE)

# C1 cluster results
c1_mrkr_res_stroma <- topTable(fit_stroma,
                               coef = "c1_vs_other",
                               n = Inf,
                               p = 1)

# C2 cluster results
c2_mrkr_res_stroma <- topTable(fit_stroma,
                               coef = "c2_vs_other",
                               n = Inf,
                               p = 1)
# C3 cluster results
c3_mrkr_res_stroma <- topTable(fit_stroma,
                               coef = "c3_vs_other",
                               n = Inf,
                               p = 1)

# C4 cluster results
c4_mrkr_res_stroma <- topTable(fit_stroma,
                               coef = "c4_vs_other",
                               n = Inf,
                               p = 1)

# List
mrkr_res_stroma <- list(
  `C1_A cluster` = c1_mrkr_res_stroma |> arrange(desc(logFC)),
  `C2_A cluster` = c2_mrkr_res_stroma |> arrange(desc(logFC)),
  `C3_A cluster` = c3_mrkr_res_stroma |> arrange(desc(logFC)),
  `C4_A cluster` = c4_mrkr_res_stroma |> arrange(desc(logFC))
)

# Supp File with DGEA results
write.xlsx(mrkr_res_stroma, rowNames = TRUE, file = "./Results/Tables/OSIRESP_stromal_subtypes_markers_stroma.xlsx")

# 2. Tumor compartment ----

# 9. Marker genes (using limma) ----

# Design matrix
design_tumor <- model.matrix( ~ 0 + cluster, data = leid_res_tumor)

fit_tumor <- lmFit(expr_tumor, design_tumor, robust = TRUE)

# Contrasts: each cluster vs the others
contrasts_mrkr_tumor <- makeContrasts(
  c1_vs_other = clusterC1 - (clusterC2 + clusterC3 + clusterC4) / 3,
  c2_vs_other = clusterC2 - (clusterC1 + clusterC3 + clusterC4) / 3,
  c3_vs_other = clusterC3 - (clusterC1 + clusterC2 + clusterC4) / 3,
  c4_vs_other = clusterC4 - (clusterC1 + clusterC2 + clusterC3) / 3,
  levels = design_tumor
)

fit_tumor <- contrasts.fit(fit_tumor, contrasts_mrkr_tumor)
fit_tumor <- eBayes(fit_tumor, trend = TRUE, robust = TRUE)

# C1 cluster results
c1_mrkr_res_tumor <- topTable(fit_tumor,
                              coef = "c1_vs_other",
                              n = Inf,
                              p = 1)

# C2 cluster results
c2_mrkr_res_tumor <- topTable(fit_tumor,
                              coef = "c2_vs_other",
                              n = Inf,
                              p = 1)
# C3 cluster results
c3_mrkr_res_tumor <- topTable(fit_tumor,
                              coef = "c3_vs_other",
                              n = Inf,
                              p = 1)

# C4 cluster results
c4_mrkr_res_tumor <- topTable(fit_tumor,
                              coef = "c4_vs_other",
                              n = Inf,
                              p = 1)

# List
mrkr_res_tumor <- list(
  `C1_A cluster` = c1_mrkr_res_tumor |> arrange(desc(logFC)),
  `C2_A cluster` = c2_mrkr_res_tumor |> arrange(desc(logFC)),
  `C3_A cluster` = c3_mrkr_res_tumor |> arrange(desc(logFC)),
  `C4_A cluster` = c4_mrkr_res_tumor |> arrange(desc(logFC))
)

# Supp File with DGEA results
write.xlsx(mrkr_res_tumor, rowNames = TRUE, file = "./Results/Tables/OSIRESP_tumorl_subtypes_markers_tumor.xlsx")

# GSVA ----

# Gene set preparation
gs_names <- cta_gs$gs_name |>
  unique()

cta_gs_list <- lapply(
  gs_names,
  FUN = function(set) {
    gene_vector <- cta_gs |>
      filter(gs_name == set &
               gene_symbol %in% colnames(expr_data)) |>
      select(gene_symbol) |>
      pull()
  }
)

# Set names to list
names(cta_gs_list) <- gs_names

# Check class
cta_gs_list |>
  class()

# Check structure
cta_gs_list |>
  str()

# Filter out gene sets with size < 15
min_len <- 15
cta_gs_list <- cta_gs_list[sapply(cta_gs_list, length) >= min_len]

# 1. Stroma compartment ----

# GSVA
gsva_stroma_par <- gsvaParam(expr_stroma, cta_gs_list)
gsva_stroma <- gsva(gsva_stroma_par, verbose = FALSE)

# 2. Tumor compartment ----

# GSVA
gsva_tumor_par <- gsvaParam(expr_tumor, cta_gs_list)
gsva_tumor <- gsva(gsva_tumor_par, verbose = FALSE)

# SAVE RESULTS ----

# Clear environment
rm(list = setdiff(
  ls(),
  c(
    "seurat_obj",
    "leid_res",
    "fit_stroma_res",
    "fit_tumor_res",
    "mrkr_res_stroma",
    "mrkr_res_tumor",
    "gsva_stroma",
    "gsva_tumor"
  )
))

save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/8_OSIRESP_stroma_unsupervised_analysis.RData")
