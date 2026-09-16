# DESCRIPTION ----
# OSIRESP cohort, including pre-treatment biopsies AND additionally collected
# post-progression biopsies.
# Summarize expression by compartment (i.e., tumor or stroma) in the multiple
# ROIs selected in the tumor sample using the median expression value per gene.

# LIBRARIES ----
library(dplyr)
library(GeomxTools)

# READ GEOMX DATA AFTER QC ----

# OSIRESP
osiresp_rb <- readRDS("./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_q3norm.rds")

# MEDIAN VALUES PER PATIENT AND COMPARMENT ----

# OSIRESP
osiresp_rb_vst <- assayDataElement(object = osiresp_rb, elt = "q3norm_vst") |>
  t()

osiresp_rb_vst <- list(pData(osiresp_rb) |>
                         select(patient_id, compartment, biopsy_type),
                       osiresp_rb_vst) |>
  bind_cols()

# Median values per gene and compartment, distinguishing between pre-treatment
# biopsies and biosies at progression
osiresp_rb_vst_median <- osiresp_rb_vst |>
  group_by(patient_id, compartment, biopsy_type) |>
  summarise_all(.funs = median)

# SAVE RDS ----

saveRDS(
  osiresp_rb_vst_median,
  "./Data/OSIRESP_RB_cohort/Processed/OSIRESP_RB_vst_medianvalues.rds"
)
