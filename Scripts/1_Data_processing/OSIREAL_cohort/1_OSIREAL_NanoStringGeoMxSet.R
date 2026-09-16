# DESCRIPTION ----
# OSIREAL cohort.
# Creation of NanoStringGeoMxSet object.

# DOCUMENTATION ----
# https://www.bioconductor.org/packages/release/bioc/html/GeomxTools.html
# https://bioconductor.riken.jp/packages/3.15/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# NOTES ----
# 1-
# The data submitted to the Gene Expression Omnibus repository include only those
# from the tumor segments included in the study after quality control
# (performed in a subsequent script).

# Since some quality control criteria apply to data from all segments, running
# the data processing scripts in this GitHub repository might yield slightly
# different results. These scripts are included in the repository solely to
# illustrate the pipeline followed for the data repository.

# 2-
# Local paths are used for the DCC files.
# DCC files can be downloaded from the Gene Omnibus Expression repository; paths
# in the script should be modified accordingly.

# 3-
# .pkc file should always be downloaded from the Gene Omnibus Expression repository;
# path in the script should be modified accordingly

# LIBRARIES ----
library(readxl)
library(openxlsx)
library(dplyr)
library(GeomxTools)
library(stringr)

# CLINICAL ANNOTATIONS ----

clinical_data <- read_xlsx("./Data/OSIREAL_cohort/Clinical_annotations/OSIREAL_clinical_annotations.xlsx")

# DCC FILE NAMES ----

# Get DCC file names (one per batch)

# Initial DCC files (paths)
path_01 <- "./Data/OSIREAL_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_01"
path_05 <- "./Data/OSIREAL_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_05"
path_R2 <- "./Data/OSIREAL_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_R2"
path_R3 <- "./Data/OSIREAL_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_R3"

# New DCC files (from resequenciation) (paths)
path_new <- "./Data/OSIREAL_cohort/Raw/OSIREAL_new_DCC"

# Initial DCC files
dcc_01 <- dir(path = path_01)
dcc_05 <- dir(path = path_05)
dcc_R2 <- dir(path = path_R2)
dcc_R3 <- dir(path = path_R3)

names(dcc_01) <- rep("dcc_01", length(dcc_01))
names(dcc_05) <- rep("dcc_05", length(dcc_05))
names(dcc_R2) <- rep("dcc_R2", length(dcc_R2))
names(dcc_R3) <- rep("dcc_R3", length(dcc_R3))

# Concatenate
dcc <- c(dcc_01, dcc_05, dcc_R2, dcc_R3)

# Check duplicates
dcc |> duplicated() |> table(useNA = "always")

# Dataframe
dcc_df <- data.frame(dcc = dcc, batch = names(dcc))

# Visualization
dcc_df |>
  View()

rm(dcc)

# New DCC files (from resequenciation)
new_dcc <- dir(path = path_new)
names(new_dcc) <- rep("dcc_resequenced", length(new_dcc))

# Remove initial DCC files corresponding to re-sequenced samples
dcc_to_exclude <- intersect(dcc_df$dcc, new_dcc)
dcc_to_exclude |> length()

# Check initial dimensions
dcc_df |>
  nrow()

# Remove
dcc_df <- dcc_df |>
  filter(!(dcc %in% dcc_to_exclude))

# Check dimensions after removing
dcc_df |>
  nrow()

# Check no intersection between initial DCC files and DCC files from resequentiation
intersect(dcc_df$dcc, new_dcc) |>
  length()

# Bind rows

# new dcc df
new_dcc_df <- data.frame(dcc = new_dcc, batch = names(new_dcc))

dcc_df <- list(dcc_df, new_dcc_df) |>
  bind_rows()

# DCC file names including complete path (one per batch)

dcc_df$batch |>
  table(useNA = "always")

dcc_df <- dcc_df |>
  mutate(
    dir = case_when(
      batch == "dcc_01" ~ paste(path_01, dcc, sep = "/"),
      batch == "dcc_05" ~ paste(path_05, dcc, sep = "/"),
      batch == "dcc_R2" ~ paste(path_R2, dcc, sep = "/"),
      batch == "dcc_R3" ~ paste(path_R3, dcc, sep = "/"),
      batch == "dcc_resequenced" ~ paste(path_new, dcc, sep = "/")
    )
  )

# PKC file ----

# PKC file name
pkc_dir <- "./Data/OSIREAL_cohort/Raw/GeoMx_Hs_CTA_v1.0.pkc"

# LABWORKSHEET ----

# Read labworksheet file
labws <- readRDS("./Data/OSIREAL_cohort/Raw/OSIREAL_labworksheet.rds")

# Number of segments in labworksheet
labws |>
  filter(`slide name` != "No Template Control") |>
  dim()

# CHECK labworksheet AND DCC FILES ----

# Check correspondence between labworksheet and DCC files
all(labws$Sample_ID %in% str_remove_all(dcc_df$dcc, ".dcc"))

all(str_remove_all(dcc_df$dcc, ".dcc") %in% labws$Sample_ID)

# Remove DCC files that are not from OSIREAL cohort
keep_bool <- str_remove_all(dcc_df$dcc, ".dcc") %in% labws$Sample_ID
keep_bool |> table()

# Check initial dimensions
dcc_df |>
  nrow()

# Remove
dcc_df <- dcc_df |>
  filter(keep_bool)

# Check dimensions after removing
dcc_df |>
  nrow()

rm(keep_bool)

# Check order of Sample_ID in labworksheet (must be the same as the one in DCC files)
all(str_remove_all(dcc_df$dcc, ".dcc") == labws$Sample_ID)

# Set same order
rownames(labws) <- labws$Sample_ID
labws <- labws[str_remove_all(dcc_df$dcc, ".dcc"), ]

# Check order
all(str_remove_all(dcc_df$dcc, ".dcc") == labws$Sample_ID)

# PHENOTYPIC DATA ----

pheno_data <- labws |>
  left_join(
    clinical_data |>
      select(
        -date_of_first_dose_of_osimertinib,-date_of_disease_progression_to_osimertinib
      ),
    by = join_by(patient_id)
  )

# Check dimensions
dim(clinical_data)
dim(labws)
dim(pheno_data)

# Check order of sample_ID in phenodata (must be the same as the one in the DCC files)
all(str_remove_all(dcc_df$dcc, ".dcc") == pheno_data$Sample_ID)
# TRUE

# Write .xlsx file (needed for the NanoStringGeoMxSet object)
write.xlsx(labws,
           "./Data/OSIREAL_cohort/Processed/OSIREAL_labworksheet.xlsx")

# path for the phenodata .xlsx file
pheno_data_dir <- "./Data/OSIREAL_cohort/Processed/OSIREAL_phenodata.xlsx"

write.xlsx(pheno_data, pheno_data_dir)

# NanoStringGeoMxSet OBJECT ----

# It is important not to change the names of the variables "scan name" and
# "slide name". If doing so, the protocolData is not created correctly and the
# QC can not be done properly
rawdata <- readNanoStringGeoMxSet(
  dccFiles = dcc_df$dir,
  pkcFiles = pkc_dir,
  phenoDataFile = pheno_data_dir,
  phenoDataSheet = "Sheet 1"
)

# Check dimensions
dim(pData(rawdata))
dim(exprs(rawdata))
dim(sData(rawdata))

pData(rawdata) |> head() |> View() # phenotypic data
exprs(rawdata) |> head() |> View() # expression and error measurements
sData(rawdata) |> head() |> View() # sample data (phenotypic data + protocol data)

# FEATURE PREPARATION FOR DOWNSTREAM ANALYSES ----

# This has to be done after the NanoStringGeoMxSet object is created, since the
# function readNanoStringGeoMxSet coerce some features to a different class

# patient_id as factor
# ecog_ps as factor
# cns_metastases as factor
# liver_metastases as factor
# early_resistance as factor

# pData
pData(rawdata) <- pData(rawdata) |>
  mutate(
    patient_id = patient_id |> factor(),
    ecog_ps_at_baseline = ecog_ps_at_baseline |>
      factor(levels = 0:3) |>
      relevel(ref = "0"),
    cns_metastases_at_baseline = cns_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    liver_metastases_at_baseline = liver_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    early_resistance = early_resistance |>
      factor() |>
      relevel(ref = "No"),
    long_response = long_response |>
      factor() |>
      relevel(ref = "No")
  )

# ecog_ps: group levels 2 and 3
pData(rawdata) <- pData(rawdata) |>
  mutate(ecog_ps_at_baseline_grouped = ecog_ps_at_baseline)

pData(rawdata)$ecog_ps_at_baseline_grouped |> levels()
levels(pData(rawdata)$ecog_ps_at_baseline_grouped) <- c("0", "1", "\u2265 2", "\u2265 2")
pData(rawdata)$ecog_ps_at_baseline_grouped |> levels()

# ecog_ps_dummie: group 0-1 and 2-3
pData(rawdata) <- pData(rawdata) |>
  mutate(ecog_ps_dummie = ecog_ps_at_baseline)
levels(pData(rawdata)$ecog_ps_dummie) <- c("< 2", "< 2", "\u2265 2", "\u2265 2")
pData(rawdata)$ecog_ps_dummie |> levels()

# SAVE RDS ----

# Labworksheet
saveRDS(labws, file = "./Data/OSIREAL_cohort/Processed/OSIREAL_labworksheet.rds")

# Raw data
saveRDS(rawdata, file = "./Data/OSIREAL_cohort/Processed/OSIREAL_NSGeoMxSet_rawdata.rds")
