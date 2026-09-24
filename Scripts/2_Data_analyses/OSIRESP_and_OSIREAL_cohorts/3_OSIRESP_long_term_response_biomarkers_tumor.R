# DESCRIPTION ----
# OSIRESP cohort, including only pre-treatment samples.
# ENET logistic regression to identify compartment-specific biomarkers associated
# with osimertinib resistance. Tumor compartment.

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

# READ GEOMX DATA AFTER QC ----

# Median expression values per tumor and compartment
data_median <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_vst_medianvalues.rds")

# Clinical annotations
pData <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_q3norm_pData.rds")

# FUNCTIONS ----
source("./Scripts/2_Data_analyses/OSIRESP_and_OSIREAL_cohorts/ENET_logistic_functions.R")


# GLOBAL PARAMETERS ----

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
unpen_vars <- "MET"
unpen_penalties <- 0

# DATA PREPARATION ----

# Tumor
model_data <- data_median |>
  ungroup() |>
  filter(compartment == "Tumor") |>
  select(-compartment) |>
  inner_join(
    pData |>
      select(patient_id, long_response, ecog_ps_dummie, egfr_mutation_type) |>
      unique()
  ) |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = ifelse(long_response == "Yes", 1, 0))

model_data$long_response |> table(useNA = "always")

# X matrix. Remove patient ID and outcome
x_matrix <- makeX(
  model_data |>
    select(
      -patient_id,-long_response,-ecog_ps_dummie,-egfr_mutation_type
    )
)

# Check dimensions
x_matrix |> dim()

# Y vector
y_vector <- model_data$long_response

# Check length
y_vector |> length()

# MODEL ----

# Set seed for replicability
set.seed(seed)

model_results <- sis_enet(
  x_matrix = x_matrix,
  y_vector = y_vector,
  unpen_vars = unpen_vars,
  unpen_penalties = unpen_penalties
)

pred_results <- predictions(
  model = model_results$model,
  lambda_opt = model_results$lambda_opt,
  newx = x_matrix
)

# ROC AUC
auc_results <- auc_estimates(model_data = model_data,
                             y_vector = y_vector,
                             pred = pred_results)

# REPEATED K-FOLD CROSS-VALIDATION ----

repcv_results <- repcv(
  model_data = model_data,
  x_matrix = x_matrix,
  y_vector = y_vector,
  unpen_vars = unpen_vars,
  unpen_penalties = unpen_penalties
)

# SAVE RESULTS ----

save.image(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/3_OSIRESP_long_term_response_biomarkers_tumor.RData")
