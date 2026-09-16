# DESCRIPTION
# Clinical data description. OSIRESP and OSIREAL cohort

# LIBRARIES ----
library(dplyr)
library(stringr)
library(gtsummary)
library(gt)

# TABLE 1 ----
# OSIRESP cohort

# OSIRESP
osiresp_qc <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_QC.rds")

# Metadata
osiresp <- pData(osiresp_qc)

# Select clinical variables
clinical_data <- osiresp |>
  mutate(
    egfr_mutation_type = egfr_mutation_type |>
      factor(),
    egfr_mutation_type_grouped = egfr_mutation_type,
    tobacco_history = tobacco_history |>
      factor() |>
      relevel(ref = "Never smoker"),
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
      factor(
        levels = c("Yes", "No", "Non-evaluable"),
        ordered = TRUE
      ),
    ecog_ps_at_baseline_grouped = ecog_ps_at_baseline,
    ecog_ps_dummie = ecog_ps_at_baseline,
    tobacco_history_grouped = tobacco_history
  ) |>
  select(
    patient_id,
    age,
    sex,
    egfr_mutation_type,
    tobacco_history,
    stage_at_baseline,
    cns_metastases_at_baseline,
    liver_metastases_at_baseline,
    ecog_ps_at_baseline,
    disease_progression_to_osimertinib,
    long_response
  ) |>
  unique()

# Table 1
clinical_data |>
  mutate(stage_at_baseline = case_when(
    str_detect(stage_at_baseline, "IV") ~ "IV",
    str_detect(stage_at_baseline, "III") ~ "III",
    str_detect(stage_at_baseline, "II") ~ "II"
  )) |>
  dplyr::rename(
    "Age" = age,
    "Sex" = sex,
    "EGFR mutation type" = egfr_mutation_type,
    "Tobacco status" = tobacco_history,
    "Stage" = stage_at_baseline,
    "CNS metastases" = cns_metastases_at_baseline,
    "Liver metastases" = liver_metastases_at_baseline,
    "ECOG PS" = ecog_ps_at_baseline,
    "Disease progression to osimertinib" = disease_progression_to_osimertinib,
    "Long-term osimertinib response" = long_response
  ) |>
  tbl_summary(
    include = c(
      "EGFR mutation type",
      Age,
      Sex,
      "Tobacco status",
      Stage,
      "CNS metastases",
      "Liver metastases",
      "ECOG PS",
      "Disease progression to osimertinib",
      "Long-term osimertinib response"
    )
  ) |>
  as_gt() |>
  gtsave(filename = "./Results/Tables/Table_1.docx")

rm(list = ls())

# TABLE S1 -----
# OSIREAL cohort

osireal_qc <- readRDS("./Data/OSIREAL_cohort/Processed/OSIREAL_QC.rds")

# Metadata
osireal <- pData(osireal_qc)

# Select clinical variables
clinical_data <- osireal |>
  mutate(
    egfr_mutation_type = egfr_mutation_type |>
      factor(),
    egfr_mutation_type_grouped = egfr_mutation_type,
    tobacco_history = tobacco_history |>
      factor() |>
      relevel(ref = "Never smoker"),
    ecog_ps_at_baseline = ecog_ps_at_baseline |>
      factor(levels = 0:3) |>
      relevel(ref = "0"),
    cns_metastases_at_baseline = cns_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    liver_metastases_at_baseline = liver_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    long_response = long_response |>
      factor(
        levels = c("Yes", "No", "Non-evaluable"),
        ordered = TRUE
      ),
    ecog_ps_at_baseline_grouped = ecog_ps_at_baseline,
    ecog_ps_dummie = ecog_ps_at_baseline,
    tobacco_history_grouped = tobacco_history
  ) |>
  select(
    patient_id,
    age,
    sex,
    egfr_mutation_type,
    tobacco_history,
    stage_at_baseline,
    cns_metastases_at_baseline,
    liver_metastases_at_baseline,
    ecog_ps_at_baseline,
    disease_progression_to_osimertinib,
    long_response
  ) |>
  unique()

# Table
clinical_data |>
  mutate(stage_at_baseline = case_when(
    str_detect(stage_at_baseline, "IV") ~ "IV",
    str_detect(stage_at_baseline, "III") ~ "III",
    str_detect(stage_at_baseline, "II") ~ "II"
  )) |>
  dplyr::rename(
    "Age" = age,
    "Sex" = sex,
    "EGFR mutation type" = egfr_mutation_type,
    "Tobacco status" = tobacco_history,
    "Stage" = stage_at_baseline,
    "CNS metastases" = cns_metastases_at_baseline,
    "Liver metastases" = liver_metastases_at_baseline,
    "ECOG PS" = ecog_ps_at_baseline,
    "Disease progression to osimertinib" = disease_progression_to_osimertinib,
    "Long-term osimertinib response" = long_response
  ) |>
  tbl_summary(
    include = c(
      "EGFR mutation type",
      Age,
      Sex,
      "Tobacco status",
      Stage,
      "CNS metastases",
      "Liver metastases",
      "ECOG PS",
      "Disease progression to osimertinib",
      "Long-term osimertinib response"
    )
  ) |>
  as_gt() |>
  gtsave(filename = "./Results/Tables/Table_S1.docx")

rm(list = ls())

# TABLE S2 ----
# OSIRESP cohort: patients with paired pre-osimertinib and post-progression
# biopsies
osiresp_rb_qc <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_QC.rds")

# Metadata
osiresp_rb <- pData(osiresp_rb_qc)

# Expression data + biopsy type annotation (pre-treatment or at progression)
expr_data <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_vst_medianvalues.rds")

# Order levels of biopsy_type
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

# Check
all(n_biopsies_stroma$patient_id %in% n_biopsies_tumor$patient_id)

# Select clinical variables
clinical_data <- osiresp_rb |>
  mutate(
    egfr_mutation_type = egfr_mutation_type |>
      factor(),
    egfr_mutation_type_grouped = egfr_mutation_type,
    tobacco_history = tobacco_history |>
      factor() |>
      relevel(ref = "Never smoker"),
    ecog_ps_at_baseline = ecog_ps_at_baseline |>
      factor(levels = 0:3) |>
      relevel(ref = "0"),
    cns_metastases_at_baseline = cns_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    liver_metastases_at_baseline = liver_metastases_at_baseline |>
      factor() |>
      relevel(ref = "No"),
    ecog_ps_at_baseline_grouped = ecog_ps_at_baseline,
    ecog_ps_dummie = ecog_ps_at_baseline,
    tobacco_history_grouped = tobacco_history,
    pfs_time_months = pfs_time_days / (365.25 / 12)
  ) |>
  select(
    patient_id,
    age,
    sex,
    egfr_mutation_type,
    tobacco_history,
    stage_at_baseline,
    cns_metastases_at_baseline,
    liver_metastases_at_baseline,
    ecog_ps_at_baseline,
    disease_progression_to_osimertinib,
    pfs_time_months
  ) |>
  unique() |>
  filter(patient_id %in% n_biopsies_tumor$patient_id) |>
  mutate(
    paired_tumor_samples = "Yes",
    paired_stromal_samples = ifelse(patient_id %in% n_biopsies_stroma$patient_id, "Yes", "No"),
    pfs_time_months = round(pfs_time_months, digits = 2)
  )

# Table
clinical_data |>
  mutate(stage_at_baseline = case_when(
    str_detect(stage_at_baseline, "IV") ~ "IV",
    str_detect(stage_at_baseline, "III") ~ "III",
    str_detect(stage_at_baseline, "II") ~ "II"
  )) |>
  dplyr::rename(
    "Patient ID" = patient_id,
    "Paired samples (tumor)" = paired_tumor_samples,
    "Paired samples (stroma)" = paired_stromal_samples,
    "Age" = age,
    "Sex" = sex,
    "EGFR mutation type" = egfr_mutation_type,
    "Tobacco status" = tobacco_history,
    "Stage" = stage_at_baseline,
    "CNS metastases" = cns_metastases_at_baseline,
    "Liver metastases" = liver_metastases_at_baseline,
    "ECOG PS" = ecog_ps_at_baseline,
    "Disease progression to osimertinib" = disease_progression_to_osimertinib,
    "Progression-free survival (months)" = pfs_time_months
  ) |>
  gt() |>
  gtsave(filename = "./Results/Tables/Table_S2.docx")

rm(list = ls())

# TABLE S3 ----
# H12O_TMA004 cohort

h12otma004_qc <- readRDS("./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_QC.rds")

# Metadata
h12otma004_qc <- pData(h12otma004_qc)

# Select clinical variables
clinical_data <- h12otma004_qc |>
  mutate(
    egfr_mutation_type = egfr_mutation_type |>
      factor(),
    tobacco_history = tobacco_history |>
      factor() |>
      relevel(ref = "Never smoker"),
    relapse = relapse |>
      factor(levels = "Yes", "No"),
    ordered = TRUE
  ) |>
  select(
    patient_id,
    age,
    sex,
    egfr_mutation_type,
    tobacco_history,
    stage_at_baseline,
    relapse
  ) |>
  unique()

# Table
clinical_data |>
  mutate(
    egfr_mutation_type = case_when(!(
      egfr_mutation_type %in% c("Exon 19 deletion", "Exon 21 L858R")
    ) ~ "Other", .default = egfr_mutation_type),
    stage_at_baseline = case_when(str_detect(stage_at_baseline, "IA") ~ "IA", .default = stage_at_baseline)
  ) |>
  dplyr::rename(
    "Age" = age,
    "Sex" = sex,
    "EGFR mutation type" = egfr_mutation_type,
    "Tobacco status" = tobacco_history,
    "Stage" = stage_at_baseline,
    "Relapse" = relapse
  ) |>
  tbl_summary(include = c("EGFR mutation type", Age, Sex, "Tobacco status", Stage, "Relapse")) |>
  as_gt() |>
  gtsave(filename = "./Results/Tables/Table_S3.docx")
