# DESCRIPTION ----
# OSIREAL cohort
# RNA data normalization: upper-quartile method.
# Varianze-stabilizing transformation.

# DOCUMENTATION ----
# https://pubmed.ncbi.nlm.nih.gov/32789507/
# https://github.com/bhattacharya-a-bt/CBCS_normalization/blob/master/CBCS_normalization_tutorial.pdf
# (RUVSeq correction is not performed)

# LIBRARIES
library(dplyr)
library(GeomxTools)
library(EDASeq)
library(DESeq2)
library(limma)
library(matrixStats)

# LOAD DATA ----

# Load target data after QC
target_data <- readRDS("./Data/OSIREAL_cohort/Processed/OSIREAL_QC.rds")
dim(target_data)

# UPPER QUANTILE NORMALIZATION AND VST TRANSFORMATION ----

# Use DESeq2 formulation to integrate raw expression
set <- newSeqExpressionSet(counts = exprs(target_data) |> round(),
                           phenoData = pData(target_data))

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

# set@assayData$normalizedCounts to access normalized counts
# assay(vsd) to access normalized and log-transformed (vst-transformed) counts

# ADD NORMALIZATION TO OBJECT ----

assayDataElement(object = target_data, elt = "q3norm") <- set@assayData$normalizedCounts
assayDataElement(object = target_data, elt = "q3norm_vst") <- assay(vsd)

# SAVE RDS ----

saveRDS(target_data,
        "./Data/OSIREAL_cohort/Processed/OSIREAL_q3norm.rds")
