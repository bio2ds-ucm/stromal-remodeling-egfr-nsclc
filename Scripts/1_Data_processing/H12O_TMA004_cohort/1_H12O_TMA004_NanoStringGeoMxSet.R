# DESCRIPTION ----
# H12O_TMA004 cohort.
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

clinical_data <- read_xlsx(
  "./Data/H12O_TMA004_cohort/Clinical_annotations/H12O_TMA004_clinical_annotations.xlsx"
)

# DCC FILE NAMES ----

# Get DCC file names (one per batch)
dcc_05 <- dir(path = "./Data/H12O_TMA004_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_05/")
dcc_R3 <- dir(path = "./Data/H12O_TMA004_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_R3")

# Concatenate
dcc <- c(dcc_05, dcc_R3)

# DCC file names including complete path (one per batch)
dcc_05_dir <- dir(path = "./Data/H12O_TMA004_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_05/", full.names = TRUE)

dcc_R3_dir <- dir(path = "./Data/H12O_TMA004_cohort/Raw/GeoMx_NGS_Pipeline_DCC_OSIRESP24_R3", full.names = TRUE)

# Concatenate
dcc_dir <- c(dcc_05_dir, dcc_R3_dir)

# PKC file ----

# PKC file name
pkc_dir <- "./Data/H12O_TMA004_cohort/Raw/GeoMx_Hs_CTA_v1.0.pkc"

# LABWORKSHEET ----

# Read labworksheet file
labws <- readRDS("./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_labworksheet.rds")

# Check dimensions
labws |>
  filter(`slide name` != "No Template Control") |>
  dim()

# CHECK LABWORKSHEET AND DCC FILES ----

all(labws$Sample_ID %in% str_remove_all(dcc, ".dcc"))

setdiff(labws$Sample_ID, str_remove_all(dcc, ".dcc"))

all(str_remove_all(dcc, ".dcc") %in% labws$Sample_ID)

setdiff(str_remove_all(dcc, ".dcc"), labws$Sample_ID) |> length()

# Remove DCC files that are not from the H12O_TMA004 cohort
keep_bool <- str_remove_all(dcc, ".dcc") %in% labws$Sample_ID
dcc <- dcc[keep_bool]
dcc_dir <- dcc_dir[keep_bool]

# Check after removing DCC files
dcc_tibble <- tibble(dcc_file = dcc, dcc_path = dcc_dir)

dcc_tibble <- dcc_tibble |>
  mutate(match = c(FALSE, TRUE)[1 + str_detect(dcc_path, as.character(dcc_file))])
dcc_tibble$match |> table(useNA = "always")

all(str_remove_all(dcc, ".dcc") %in% labws$Sample_ID)

rm(keep_bool)

# Check order of Sample_ID in labworksheet (must be the same as the one in DCC files)
all(str_remove_all(dcc, ".dcc") == labws$Sample_ID)

# Set same order
rownames(labws) <- labws$Sample_ID
labws <- labws[str_remove_all(dcc, ".dcc"), ]

# Check order
all(str_remove_all(dcc, ".dcc") == labws$Sample_ID)

# PHENOTYPIC DATA ----

pheno_data <- labws |>
  left_join(
    clinical_data |>
      select(-date_of_surgery, -date_of_last_follow_up),
    by = join_by(patient_id)
  )

# Check dimensions
dim(clinical_data)
dim(labws)
dim(pheno_data)

# Check order of sample_ID in phenoData (must be the same as the one in the DCC files)
all(str_remove_all(dcc, ".dcc") == pheno_data$Sample_ID)

# Write .xlsx file (needed for the NanoStringGeoMxSet object)
write.xlsx(labws,
           "./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_labworksheet.xlsx")

# path for the phenoData .xlsx file
pheno_data_dir <- "./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_phenoData.xlsx"

write.xlsx(pheno_data, pheno_data_dir)

# NanoStringGeoMxSet OBJECT ----

# It is important not to change the names of the variables "scan name" and
# "slide name". If doing so, the protocolData is not created correctly and the
# QC can not be done properly
rawdata <- readNanoStringGeoMxSet(
  dccFiles = dcc_dir,
  pkcFiles = pkc_dir,
  phenoDataFile = pheno_data_dir,
  phenoDataSheet = "Sheet 1"
)

# Check dimensions
dim(pData(rawdata))
dim(exprs(rawdata))
dim(sData(rawdata))
# # 457 segments

pData(rawdata) |> head() |> View() # phenotypic data
exprs(rawdata) |> head() |> View() # expression and error measurements
sData(rawdata) |> head() |> View() # sample data (phenotypic data + protocol data)

# SAVE RDS ----

# Labworksheet
saveRDS(labws, file = "./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_labworksheet.rds")

# Raw data
saveRDS(rawdata, file = "./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_NSGeoMxSet_rawdata.rds")
