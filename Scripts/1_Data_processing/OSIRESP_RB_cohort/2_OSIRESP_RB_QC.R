# DESCRIPTION ----
# OSIRESP cohort, including pre-treatment biopsies AND additionally collected
# post-progression biopsies.
# Quality control of tumor samples following Bruker NanoString recommendations

# DOCUMENTATION ----
# https://bioconductor.riken.jp/packages/3.15/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html


# NOTES ----
# 1-
# The data submitted to the Gene Expression Omnibus repository include only those
# from the tumor segments included in the study quality control.

# Since some quality control criteria apply to data from all segments, running
# the data processing scripts in this GitHub repository might yield slightly
# different results. These scripts are included in the repository solely to
# illustrate the pipeline followed for the data repository.


# 2-
# Quality control is performed in several steps:
#
# 1. Segment QC, remove ROIs or segments that do not fulfill criteria
# 2. Probe QC, remove probes that do not fulfill criteria
# 3. Target aggregation
# 4. Limit of Quantification (LOQ) per segment
#       discard segments with levels of expression consistently below LOQ

# LIBRARIES ----
library(dplyr)
library(GeomxTools)

# DATA PREPARATION ----

# Load NanoStringGeoMxSet object
rawdata <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_NSGeoMxSet_rawdata.rds")

# Shift all counts in expression matrix by 1 (shift all, not only those that are 0)
data <- shiftCountsOne(rawdata, useDALogic = FALSE)

# 1. SEGMENT QC ----

# QC parameters ----

qc_params <- list(
  minSegmentReads = 1000,
  # Minimum number of reads
  percentTrimmed = 80,
  # Minimum % of reads trimmed
  percentStitched = 80,
  # Minimum % of reads stitched
  percentAligned = 80,
  # Minimum % of reads aligned
  percentSaturation = 30,
  # Minimum sequencing saturation
  minNegativeCount = 5,
  # Minimum negative control counts
  maxNTCCount = 10000,
  # Maximum counts observed in NTC well
  minNuclei = 100,
  # Minimum # of nuclei estimated
  minArea = 10000
)        # Minimum segment area

# Flag segments ----
data <- setSegmentQCFlags(data, qcCutoffs = qc_params)

# Save the results of QC
protocolData(data)[["QCFlags"]] |> View()
qc_results <- protocolData(data)[["QCFlags"]]

flag_columns <- colnames(qc_results)

# Flags
qc_summary <- data.frame(Pass = colSums(!qc_results[, flag_columns]),
                         Warning = colSums(qc_results[, flag_columns]))
qc_summary

# Total number of samples flagged
qc_results$QCStatus <- apply(qc_results, 1L, function(x) {
  ifelse(sum(x) == 0L, "Pass", "Warning")
})

qc_summary["Total flags", ] <-
  c(sum(qc_results[, "QCStatus"] == "Pass"), sum(qc_results[, "QCStatus"] == "Warning"))

qc_summary

# Negative negometric mean ----

# Calculate the negative geometric means for each module
negativeGeoMeans <- esBy(
  negativeControlSubset(data),
  GROUP = "Module",
  FUN = function(x) {
    assayDataApply(x,
                   MARGIN = 2,
                   FUN = ngeoMean,
                   elt = "exprs")
  }
)

protocolData(data)[["NegGeoMean"]] <- negativeGeoMeans

# Explicitly copy the Negative geoMeans from sData to pData
pkcs <- annotation(rawdata)
modules <- gsub(".pkc", "", pkcs)
negCols <- paste0("NegGeoMean_", modules)
pData(data)[, negCols] <- sData(data)[["NegGeoMean"]]

# Detatch neg_geomean columns ahead of aggregateCounts call
pData(data) <- pData(data)[, !colnames(pData(data)) %in% negCols]

# Discard samples not passing the QC ----

# Dimensions before segment QC
dim(data)

data <- data[, qc_results$QCStatus == "Pass"]

# Dimensions after segment QC
dim(data)

# Samples not passing the filters
qc_results[which(qc_results$QCStatus != "Pass"), ] |>
  View()

# 2. PROBE QC ----

data <- setBioProbeQCFlags(
  data,
  qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20),
  removeLocalOutliers = TRUE
)

probeqc_results <- fData(data)[["QCFlags"]]

# Define QC table for probe QC
qc_df <- data.frame(
  Passed = sum(rowSums(probeqc_results[, -1]) == 0),
  Global = sum(probeqc_results$GlobalGrubbsOutlier),
  Local = sum(
    rowSums(probeqc_results[, -2:-1]) > 0 &
      !probeqc_results$GlobalGrubbsOutlier
  )
)

# Visualise how many probes do not pass the filter
qc_df

# Exclude outlier probes
probeQCpassed <- subset(data, fData(data)[["QCFlags"]][, c("LowProbeRatio")] == FALSE &
                          fData(data)[["QCFlags"]][, c("GlobalGrubbsOutlier")] == FALSE)

# Dimensions before probe QC
dim(data)

# Dimensions after probe QC
dim(probeQCpassed)

# Remove flagged probes
data <- probeQCpassed

# 3. TARGET AGGREGATION ----

target_data <- aggregateCounts(data)

# Number of targets
length(unique(featureData(data)[["TargetName"]]))
dim(target_data)

# Visualization of expression matrix after target aggregation
# (genes in rows, samples in columns)
exprs(target_data)[1:5, 1:2]

# Manually include the negative control probe, for downstream use
negativeprobefData <- subset(fData(target_data), CodeClass == "Negative")
neg_probes <- unique(negativeprobefData$TargetName)

saveRDS(
  neg_probes,
  "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_QC_negative_probes.rds"
)

# 4. LIMIT OF QUANTIFICATION AND SEGMENT FILTERING ----

# 4.1. Limit of quantification (LOQ) ----

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_data))
for (module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if (all(vars[1:2] %in% colnames(pData(target_data)))) {
    LOQ[, module] <-
      pmax(minLOQ, pData(target_data)[, vars[1]] *
             pData(target_data)[, vars[2]] ^ cutoff)
  }
}

pData(target_data)$LOQ <- LOQ

# Filter out  segments with abnormally low signal
# Determine the number of genes detected in each segment across the dataset
LOQ_mat <- c()
for (module in modules) {
  ind <- fData(target_data)$Module == module
  Mat_i <- t(esApply(
    target_data[ind, ],
    MARGIN = 1,
    FUN = function(x) {
      x > LOQ[, module]
    }
  ))
  LOQ_mat <- rbind(LOQ_mat, Mat_i)
}

# Ensure ordering since this is stored outside of the geomxSet
LOQ_mat <- LOQ_mat[fData(target_data)$TargetName, ]

# Cheeck dimensions
dim(target_data)
dim(LOQ_mat)

# Check order
all.equal(rownames(pData(target_data)), colnames(LOQ_mat))

# 4.2. Segment filtering ----

# Save detection rate information to pheno data
# Number of genes passing the LOQ
pData(target_data)$GenesDetected <- colSums(LOQ_mat, na.rm = TRUE)

# Proportion of genes passing the LOQ
pData(target_data)$GeneDetectionRate <- pData(target_data)$GenesDetected / nrow(target_data)

# Filter out segments with < 10% of genes above LOQ

# Dimensions before segment filtering
dim(target_data)

# Segment filtering
target_data <- target_data[, pData(target_data)$GeneDetectionRate >= 0.10]

# Dimensions after segment filtering
dim(target_data)

# SAVE RESULTS ----
saveRDS(target_data,
        "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_QC.rds")
