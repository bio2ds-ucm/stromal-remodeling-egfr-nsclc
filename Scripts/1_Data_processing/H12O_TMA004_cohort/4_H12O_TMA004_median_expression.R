# DESCRIPTION ----
# H12O_TMA004 cohort.
# Summarize expression by compartment (i.e., tumor or stroma) in the multiple
# ROIs selected within a tumor using the median expression value per gene.

# LIBRARIES ----
library(dplyr)
library(GeomxTools)

# READ GEOMX DATA AFTER QC ----

expr_data <- readRDS("./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_q3norm.rds")

# MEDIAN VALUES PER PATIENT AND COMPARMENT ----

expr_data_vst <- assayDataElement(object = expr_data, elt = "q3norm_vst") |>
  t()

expr_data_vst <- list(pData(expr_data) |>
                        select(patient_id, compartment),
                      expr_data_vst) |>
  bind_cols()

expr_data_vst_median <- expr_data_vst |>
  group_by(patient_id, compartment) |>
  summarise_all(.funs = median)

expr_data_vst_median |>
  group_by(compartment) |>
  summarise(n = n())

expr_data_vst |>
  group_by(patient_id) |>
  summarise(n = n()) |>
  View()

# SAVE RDS ----

saveRDS(
  expr_data_vst_median,
  "./Data/H12O_TMA004_cohort/Processed/H12O_TMA004_vst_medianvalues.rds"
)