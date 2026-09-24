# DESCRIPTION ----
# OSIRESP cohort, including pre-treatment biopsies AND additionally collected
# post-progression biopsies.
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

# 4-
# The labworksheet used here does include post-progression biopsies taken from a
# subset of patients (unlike the one used in the processing of OSIRESP cohort data
# for the study’s primary analyses, which consider only baseline or pre-treatment biopsies).

# LIBRARIES -----
library(readxl)
library(openxlsx)
library(dplyr)
library(GeomxTools)
library(stringr)

# CLINICAL ANNOTATIONS ----

clinical_data <- read_xlsx(
  "./Data/OSIRESP_RB_cohort/Clinical_annotations/OSIRESP_RB_clinical_annotations.xlsx"
)

# DCC FILE NAMES ----

# Get DCC file names (one per batch)

# OSIRESP_I
dcc_I_01 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220202_Osiresp1_Placas1-2-3")
dcc_I_02 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220117_Osiresp2_Placas4-5-6")
dcc_I_03 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220120_Osiresp3_Placas7-8-9")
dcc_I_04 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220202_Osiresp4_Placas10-11-12")

# OSIRESP II
dcc_II_01 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP II/DCC-20230516_OsirespV_Placa1")
dcc_II_02 <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP II/DCC-20230524_OsirespV_Placas2-3")

# Concatenate
dcc <- c(dcc_I_01, dcc_I_02, dcc_I_03, dcc_I_04, dcc_II_01, dcc_II_02)

# DCC file names including complete path (one per batch)

# OSIRESP_I
dcc_I_01_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220202_Osiresp1_Placas1-2-3", full.names = TRUE)
dcc_I_02_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220117_Osiresp2_Placas4-5-6", full.names = TRUE)
dcc_I_03_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220120_Osiresp3_Placas7-8-9", full.names = TRUE)
dcc_I_04_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP I/DCC-20220202_Osiresp4_Placas10-11-12", full.names = TRUE)


# OSIRESP II
dcc_II_01_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP II/DCC-20230516_OsirespV_Placa1", full.names = TRUE)
dcc_II_02_dir <- dir(path = "./Data/OSIRESP_cohort/Raw/OSIRESP II/DCC-20230524_OsirespV_Placas2-3", full.names = TRUE)

# Concatenate
dcc_dir <- c(
  dcc_I_01_dir,
  dcc_I_02_dir,
  dcc_I_03_dir,
  dcc_I_04_dir,
  dcc_II_01_dir,
  dcc_II_02_dir
)

# PKC FILE ----

# PKC file name
pkc_dir <- "./Data/OSIRESP_cohort/Raw/GeoMx_Hs_CTA_v1.0.pkc"

# LABWORKSHEET ----

# Read labworksheet file
labws <- readRDS("./Data/OSIRESP_RB_cohort/Raw/OSIRESP_RB_labworksheet.rds")

# Check dimensions
labws |>
  filter(`slide name` != "No Template Control") |>
  dim()

# Check correspondence between labworksheet and DCC files
all(labws$Sample_ID %in% str_remove_all(dcc, ".dcc"))
setdiff(labws$Sample_ID, str_remove_all(dcc, ".dcc"))

all(str_remove_all(dcc, ".dcc") %in% labws$Sample_ID)
setdiff(str_remove_all(dcc, ".dcc"), labws$Sample_ID) |> length()

# Remove DCC files that are not from OSIRESP cohort
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
      select(
        -date_of_first_dose_of_osimertinib,-date_of_disease_progression_to_osimertinib
      ),
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
           "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_labworksheet.xlsx")

# Path for the phenoData .xlsx file
pheno_data_dir <- "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_phenoData.xlsx"

# Write .xlsx file with phenoData
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
# # 1324       segments

pData(rawdata) |> head() |> View() # phenotypic data
exprs(rawdata) |> head() |> View() # expression and error measurements
sData(rawdata) |> head() |> View() # sample data (phenotypic data + protocol data)

# Add pre-treatment or at progression biopsy to pData
pData(rawdata) |> class()
pData(rawdata) <- pData(rawdata) |>
  mutate(biopsy_type = case_when(str_detect(`slide name`, "RB_") ~ "at progression", .default = "pre-treatment"))

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
saveRDS(labws, file = "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_labworksheet.rds")

# Raw data
saveRDS(rawdata, file = "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_NSGeoMxSet_rawdata.rds")
