# DESCRIPTION ----
# OSIREAL cohort.
# Obtaining raw pseudo-bulk counts per tumor, followed by normalization and
# variance-stabilizing transformation.

# LIBRARIES ----
library(dplyr)
library(GeomxTools)
library(EDASeq)
library(DESeq2)
library(limma)
library(matrixStats)

# LOAD RAW EXPRESSION DATA ----

# Load target data after QC
target_data <- readRDS("./Data/OSIRESP_cohort/Processed/OSIREAL_QC.rds")
dim(target_data)

expr_data <- exprs(target_data)
pData <- pData(target_data)

# RAW PSEUDO-BULK COUNTS ----

# Transpose
expr_data_t <- expr_data |>
  t()

# Add patient ID
expr_data_t <- pData |>
  select(patient_id, compartment) |>
  bind_cols(expr_data_t)

# Add up the raw counts per tumor (that is, from all ROIs, without distinguishing
# by compartment)
expr_pseudobulk <- expr_data_t |>
  select(-compartment) |>
  group_by(patient_id) |>
  summarise_all(.funs = sum)

# Check dimensions
dim(expr_pseudobulk)

pData_pat <- pData |>
  select(patient_id, long_response) |>
  unique()

# Check dimensions
dim(pData_pat)

# NORMALIZATION AND VST TRANSFORMATION ----

# Check order
all(expr_pseudobulk$patient_id == pData_pat$patient_id)
# FALSE

# Set the same order
expr_pseudobulk <- expr_pseudobulk |>
  mutate(patient_id = factor(patient_id, levels = pData_pat$patient_id)) |>
  arrange(patient_id) |>
  mutate(patient_id = patient_id |> as.character())

# Check after reordering
all(expr_pseudobulk$patient_id == pData_pat$patient_id)

pseudobulk_counts <- expr_pseudobulk |>
  select(-patient_id) |>
  t()

# Use DESeq2 formulation to integrate raw expression
set <- newSeqExpressionSet(counts = pseudobulk_counts |> round(), phenoData = pData_pat)

# phenoData(set)@data to accesss phenoData in set

# Upper quantile normalization (Bullard 2010)
set <- betweenLaneNormalization(set, which = "upper")

# Size factor estimation
dds <- DESeqDataSetFromMatrix(counts(set), colData = pData(set), design = ~
                                1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
cts <- counts(dds, normalized = TRUE)
disp <- pmax((rowVars(cts) - rowMeans(cts)), 0) / rowMeans(cts) ^ 2
mcols(dds)$dispGeneEst <- disp
dds <- estimateDispersionsFit(dds, fitType = "mean")

# Transformation to the log space with a variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

vst_data <- assay(vsd) |>
  t()

pseudobulk_data <- pData(set) |>
  bind_cols(vst_data)

# Subset patients with both tumor and stroma RNA data
n_compart <- expr_data_t |>
  select(patient_id, compartment) |>
  unique() |>
  group_by(patient_id) |>
  summarise(n = n())

pat_subset <- n_compart |>
  filter(n == 2)

pseudobulk_data <- pseudobulk_data |>
  filter(patient_id %in% pat_subset$patient_id)

# SAVE RESULTS ----

saveRDS(pseudobulk_data,
        "./Data/OSIREAL_cohort/Processed/OSIREAL_pseudobulk.rds")
