# DESCRIPTION ----
# OSIRESP cohort, including only pre-treatment samples.
# ENET Cox regression to identify compartment-specific biomarkers associated
# with osimertinib resistance. Stromal compartment.

# NOTES ----
# 1-
# This script was run in a workstation using multiple cores (parallelization)
# due to computational demands.
# Keep in mind that it takes a long time to run.

# LIBRARIES ----
library(dplyr)
library(tidyr)
library(furrr)
library(SIS)
library(glmnet)
library(pROC)
library(survival)
library(timeROC)

# READ GEOMX DATA AFTER QC ----

# Median expression values per tumor and compartment
data_median <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_vst_medianvalues.rds")

# Clinical annotations
pData <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_q3norm_pData.rds")

# FUNCTIONS ----
source("./Scripts/2_Data_analyses/OSIRESP_and_OSIREAL_cohorts/ENET_cox_functions.R")

# GLOBAL PARAMETLRS ----

# Seed
seed <- 18501850

# Number of folds
nfolds <- 5

# Length of lambda grid
nlambda <- 100

# Alpha grid
alpha_grid <- seq(from = 0.2, to = 0.8, by = 0.2)

# Number of repetitions for hyperparameter tuning (repeated cross-validation)
cvrep_tuning <- 100

# Number of repetitions for internal validation (repeated cross-validation)
cvrep <- 100

# Number of cores for parallelization
ncores <- 10

# Unpenalized predictors
unpen_vars <- NULL
unpen_penalties <- NULL

# Times (for time dependent ROC AUC)
times <- c(10, 24)

# DATA PREPARATION ----

# Disease progression formatting
# Yes = 1, No = 0
pData$disease_progression_to_osimertinib |> table(useNA = "always")

pData <- pData |>
  mutate(
    disease_progression_to_osimertinib = case_when(disease_progression_to_osimertinib == "Yes" ~ 1, .default = 0)
  )

pData$disease_progression_to_osimertinib |> table(useNA = "always")

# Stroma
model_data <- data_median |>
  ungroup() |>
  filter(compartment == "Stroma") |>
  select(-compartment) |>
  inner_join(
    pData |>
      select(
        patient_id,
        disease_progression_to_osimertinib,
        pfs_time_months,
        ecog_ps_dummie,
        egfr_mutation_type
      ) |>
      unique()
  )

model_data |> dim()

# X matrix. Remove patient ID and outcome
x_matrix <- makeX(
  model_data |>
    select(
      -patient_id,-disease_progression_to_osimertinib,
      -pfs_time_months,-ecog_ps_dummie,-egfr_mutation_type
    )
)

# Check dimensions
x_matrix |> dim()

# Y vector
y_surv <- Surv(
  time = model_data$pfs_time_months,
  event = model_data$disease_progression_to_osimertinib
)

# Check dimensions
y_surv |> dim()

# MODEL ----

# Set seed for replicability
set.seed(seed)

model_results <- sis_enet_cox(
  x_matrix = x_matrix,
  y_surv = y_surv,
  unpen_vars = unpen_vars,
  unpen_penalties = unpen_penalties
)

# Linear predictor
pred_results <- predictions_cox(
  model = model_results$model,
  lambda_opt = model_results$lambda_opt,
  newx = x_matrix
)

# C index estimates
c_index_results <- c_index_estimates(
  model = model_results$model,
  model_data = model_data,
  lambda_opt = model_results$lambda_opt,
  x_matrix = x_matrix,
  y_surv = y_surv
)

# AUC time dependent ROC
timeroc_results <- timeroc_estimates(
  model_data = model_data,
  y_times = y_surv[, 1],
  y_events = y_surv[, 2],
  marker = pred_results,
  cause = 1,
  times = times
)

# Predicted survival probabilities at specified timepoints
surv_prob <- survfit_prob(
  model = model_results$model,
  lambda_opt = model_results$lambda_opt,
  x_matrix = x_matrix,
  y_surv = y_surv,
  newx = x_matrix,
  times = times
)

colnames(surv_prob) <- model_data$patient_id
rownames(surv_prob) <- paste0("time", times)

# REPEATED K-FOLD CROSS-VALIDATION ----

# To obtain estimates of the time-dependent ROC curve AUC (i.e., to quantify
# the ability of the tumor transcriptome to discriminate osimertinib resistance)
repcv_results <- repcv_cox(
  model_data = model_data,
  x_matrix = x_matrix,
  y_surv = y_surv,
  unpen_vars = unpen_vars,
  unpen_penalties = unpen_penalties
)

# SAVE RESULTS ----

save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/2_OSIRESP_resistance_biomarkers_stroma.RData")
