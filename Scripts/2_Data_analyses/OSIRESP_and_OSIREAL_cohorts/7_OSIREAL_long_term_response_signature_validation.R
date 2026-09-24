# DESCRIPTION ----

# Validation of both the spatially resolved gene signature and the pseudobulk
# signature for discrimination of long-term osimertinib response using the
# independent cohort of patients with advanced EGFR-mutant NSCLC OSIREAL.

# The transcriptomic data from the OSIREAL cohort used for validation were processed
# in the same way as those from OSIRESP for the development of each of the gene signatures:
#
# - Spatially resolved signature: data with compartment-specific resolution
# (i.e., tumor or stroma), median VST-transformed counts per gene and compartment
# for each patient

# - Pseudo-bulk signature: VST-transformed pseudo-bulk counts per gene for each patient

# LIBRARIES ----
library(dplyr)
library(janitor)
library(GeomxTools)
library(randomForestSRC)
library(pROC)
library(stringr)

# SPATIALLY RESOLVED GENE SIGNATURE ----

rf_final <- readRDS(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/6_OSIRESP_long_term_response_signature_spatial.rds")

# Data for validation (OSIRESP cohort) ----

# 1. Loading ----
test_data <- readRDS("./Data/OSIREAL_cohort/Processed/OSIREAL_q3norm.rds")

test_data_median <- readRDS("./Data/OSIREAL_cohort/Processed/OSIREAL_vst_medianvalues.rds")

# 2. Preparation ----
model_test_data_tumor <- test_data_median |>
  ungroup() |>
  filter(compartment == "Tumor") |>
  select(-compartment) |>
  inner_join(
    pData(test_data) |>
      select(
        patient_id,
        ecog_ps_dummie,
        egfr_mutation_type,
        disease_progression_to_osimertinib,
        pfs_time_months,
        long_response
      ) |>
      unique()
  )

# Check dimensions
model_test_data_tumor |> dim()

model_test_data_stroma <- test_data_median |>
  ungroup() |>
  filter(compartment == "Stroma") |>
  select(-compartment) |>
  inner_join(
    pData(test_data) |>
      select(
        patient_id,
        ecog_ps_dummie,
        egfr_mutation_type,
        disease_progression_to_osimertinib,
        pfs_time_months,
        long_response
      ) |>
      unique()
  )

# Check dimensions
model_test_data_stroma |> dim()

# Join
model_test_data <- model_test_data_tumor |>
  inner_join(
    model_test_data_stroma,
    by = join_by(
      patient_id,
      ecog_ps_dummie,
      egfr_mutation_type,
      long_response,
      disease_progression_to_osimertinib,
      pfs_time_months,
    ),
    suffix = c("_tumor", "_stroma")
  ) |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = long_response |> factor())

model_test_data |> dim()

# Special characters in gene names
colnames(model_test_data) <- make_clean_names(colnames(model_test_data), case = "screaming_snake")

# Recover lowercase
model_test_data <- model_test_data |>
  rename(
    disease_progression_to_osimertinib = DISEASE_PROGRESSION_TO_OSIMERTINIB,
    pfs_time_months = PFS_TIME_MONTHS,
    long_response = LONG_RESPONSE,
    patient_id = PATIENT_ID,
    ecog_ps_dummie = ECOG_PS_DUMMIE,
    egfr_mutation_type = EGFR_MUTATION_TYPE
  )

colnames(model_test_data) <- str_replace_all(colnames(model_test_data),
                                             c("TUMOR" = "tumor", "STROMA" = "stroma"))

model_test_data_predictors <- model_test_data[, rf_final$xvar.names]

# Gene signature predictions ----

pred_test <- predict.rfsrc(rf_final, newdata = model_test_data_predictors)

pred_test |> str()

advanced_pred <- tibble(
  patient_id = model_test_data$patient_id,
  disease_progression_to_osimertinib = model_test_data$disease_progression_to_osimertinib,
  pfs_time_months = model_test_data$pfs_time_months,
  long_response = model_test_data$long_response,
  pred = pred_test$predicted[, "Yes"]
) |>
  arrange(patient_id)

# Clear environment
rm(list = setdiff(ls(), c("rf_final", "advanced_pred")))

# PSEUDO-BULK GENE SIGNATURE ----

rf_final_pseudo <- readRDS(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/6_OSIRESP_long_term_response_signature_pseudobulk.rds"
)

# Data for validation (OSIRESP cohort, pseudobulk) ----

# 1. Loading ----

test_data_pseudobulk <- readRDS("./Data/OSIREAL_cohort/Processed/OSIREAL_pseudobulk.rds")

# 2. Preparation ----

test_data_pseudobulk <- test_data_pseudobulk |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = long_response |> factor())

# Special characters in gene names
colnames(test_data_pseudobulk) <- janitor::make_clean_names(colnames(test_data_pseudobulk), case = "screaming_snake")

# Recover lowercase
test_data_pseudobulk <- test_data_pseudobulk |>
  rename(long_response = LONG_RESPONSE, patient_id = PATIENT_ID)


# Gene signature predictions ----

pred_test_pseudo <- predict.rfsrc(rf_final_pseudo, newdata = test_data_pseudobulk)

pred_test_pseudo |> str()

advanced_pred_pseudo <- tibble(
  patient_id = test_data_pseudobulk$patient_id,
  long_response = test_data_pseudobulk$long_response,
  pred = pred_test_pseudo$predicted[, "Yes"]
) |>
  arrange(patient_id)

# DISCRIMINATION CAPACITY OF SPATIALLY RESOLVED AND PSEUDOBULK SIGNATURES ----

seed <- 123456
n.boot <- 10000

# ROC AUC for the spatially resolved gene signature
advanced_auc <- roc(
  response = advanced_pred$long_response,
  predictor = advanced_pred$pred,
  direction = "<"
)

# Bootstrap confidence interval
set.seed(seed)
advanced_auc_ci <- ci.auc(advanced_auc, method = "boot", boot.n = n.boot)

# ROC AUC for the pseudobulk gene signature
advanced_auc_pseudo <- roc(
  response = advanced_pred_pseudo$long_response,
  predictor = advanced_pred_pseudo$pred,
  direction = "<"
)

# Bootstrap confidence interval
set.seed(seed)
advanced_auc_ci_pseudo <- ci.auc(advanced_auc_pseudo, method = "boot", boot.n = n.boot)

# Clear environment
rm(list = setdiff(
  ls(),
  c(
    "advanced_pred",
    "early_pred",
    "advanced_auc",
    "advanced_auc_ci",
    "advanced_pred_pseudo",
    "advanced_auc_pseudo",
    "advanced_auc_ci_pseudo"
  )
))

# SAVE RESULTS ----

save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/7_OSIRESP_long_term_response_signature_validation.RData")
