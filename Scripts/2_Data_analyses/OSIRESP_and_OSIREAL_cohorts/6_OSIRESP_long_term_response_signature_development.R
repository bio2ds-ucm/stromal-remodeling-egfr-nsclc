# DESCRIPTION ----
# OSIRESP cohort, including only pre-treatment samples.
# Development of a spatially resolved gene signature with compartment resolution
# (i.e., tumor or stromal compartment) for discrimination of long-term osimertinib
# response, considering biomarkers from both the tumor and stromal compartments.

# Stringent feature selection using the variable priority algorithm, following by
# the Random forest machine learning algorithm.

# In addition, development of an analogous gene signature based on pseudobulk
# gene expression to assess the relevance of spatial (i.e., compartment-specific)
# resolution in RNA profiling.

# LIBRARIES ----
library(dplyr)
library(tidyr)
library(varPro)
library(randomForestSRC)
library(pROC)
library(stringr)
library(ggpubr)
library(janitor)

# READ GEOMX DATA AFTER QC ----

# Median expression values per tumor and compartment
data_median <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_vst_medianvalues.rds")

# Clinical annotations
pData <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_q3norm_pData.rds")

# FUNCTIONS ----
source("./Scripts/2_Data_analyses/OSIRESP_and_OSIREAL_cohorts/RF_functions.R")

# DATA PREPARATION ----

# Tumor
model_data_tumor <- data_median |>
  ungroup() |>
  filter(compartment == "Tumor") |>
  select(-compartment) |>
  inner_join(pData |>
               select(patient_id,
                      long_response,
                      pfs_time_months,
                      disease_progression_to_osimertinib,
                      ecog_ps_dummie,
                      egfr_mutation_type) |>
               unique()) |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = long_response |> factor()) |>
  as.data.frame()

model_data_tumor |> dim()

# Stroma
model_data_stroma <- data_median |>
  ungroup() |>
  filter(compartment == "Stroma") |>
  select(-compartment) |>
  inner_join(pData |>
               select(patient_id,
                      long_response,
                      pfs_time_months,
                      disease_progression_to_osimertinib,
                      ecog_ps_dummie,
                      egfr_mutation_type) |>
               unique()) |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = long_response |> factor()) |>
  as.data.frame()

model_data_stroma |> dim()

# Join
model_data <- model_data_tumor |>
  inner_join(model_data_stroma,
             by = join_by(patient_id, long_response, 
                          pfs_time_months,
                          disease_progression_to_osimertinib,
                          ecog_ps_dummie, egfr_mutation_type),
             suffix = c("_tumor", "_stroma"))

model_data |> dim()

model_data_final <- model_data |>
  select(-patient_id,
         -pfs_time_months,
         -disease_progression_to_osimertinib,
         -ecog_ps_dummie,
         -egfr_mutation_type) 

model_data_final |> dim()

# Special characters in gene names
colnames(model_data_final) <- janitor::make_clean_names(colnames(model_data_final),
                                                        case = "screaming_snake")

# Recover lowercase
model_data_final <- model_data_final |>
  dplyr::rename(long_response = LONG_RESPONSE)

colnames(model_data_final) <- str_replace_all(colnames(model_data_final), 
                                              c("TUMOR" = "tumor", "STROMA" = "stroma"))

# VARIABLE SELECTION USING VARPRO -----

# Documentation: https://cran.r-project.org/web/packages/varPro/index.html

# Seed
seed <- 123456

# Number of observations
nobs <- model_data_final |> nrow()

max_ntree <- 1000

# Feature selection with varPro
set.seed(seed)
options(rf.cores = 1)
options(mc.cores = 1)
varpro_results <- varPro:::varpro(long_response ~ .,
                                  data = model_data_final,
                                  method = "randomForestSRC",
                                  use.rfq = TRUE,
                                  nvar = nobs - 1,
                                  ntree = max_ntree,
                                  seed = seed,
                                  parallel = FALSE,
                                  cores = 1,
                                  papply = lapply)

set.seed(seed)
varpro_imp <- varPro:::importance(varpro_results,
                                  papply = lapply)

# Variable selected 
varpro_vars <- varpro_imp$unconditional |>
  rownames()

# Features to keep
tokeep_vars <- "MET_tumor"

topvars <- c(varpro_vars, tokeep_vars) |> unique()

model_data_subset <- model_data[, c("long_response", topvars)]

model_data_subset |> dim()

n_varpro_vars <- topvars |> length()

# Number of variables selected
n_varpro_vars

# GLOBAL PARAMETERS -----

# Number of features
# variables to keep in pre-filtering
n_top <- n_varpro_vars
n_top_sqrt <- n_top |> sqrt() 

nodesize_grid <- floor(c(0.02, 0.05, 0.10, 0.15)*nobs)

mtry_grid <- floor(c(n_top_sqrt/2, n_top_sqrt, n_top_sqrt*1.5, 
                     n_top_sqrt*2, n_top/3, n_top/2, 
                     n_top)) |> unique() |> sort() 

# Split rule
splitrule <- "gini"

# SPATIALLY RESOLVED PREDICTION MODEL ----

# Hyperparameter tuning ----

# 1. mtry and nodesize ----

# Set seed
# set.seed(seed)
rf_tuning_results <- rf_tuning(data = model_data_subset,
                               mtry_grid = mtry_grid,
                               nodesize_grid = nodesize_grid,
                               max_ntree = max_ntree,
                               seed = seed,
                               splitrule = splitrule)

rf_tuning_results

# 2. ntree ----

# Set seed
# set.seed(seed)
# ntree
# (using optimal values for mtry and nodesize)
rf_tuning_ntree_results <- rf_tuning_ntree(data = model_data_subset,
                                           mtry = rf_tuning_results$mtry,
                                           nodesize = rf_tuning_results$nodesize,
                                           max_ntree = max_ntree,
                                           seed = seed,
                                           splitrule = splitrule)

rf_tuning_ntree_results

plot(rf_tuning_ntree_results$rf_object)

ntree_opt <- rf_tuning_ntree_results$ntree_opt[, "ntree"] |>
  as.numeric()

# Final RF model ----

# Set seed
# set.seed(seed)
spatial_rf <- imbalanced(long_response ~ ., 
                         data = model_data_subset,
                         method = "rfq",
                         splitrule = splitrule,
                         # perf.type = "g.mean",
                         ntree = ntree_opt, 
                         mtry = rf_tuning_results$mtry,
                         nodesize = rf_tuning_results$nodesize,
                         nsplit = 10,
                         # importance = "permute",
                         seed = seed)

# OOB predictions
spatial_rf_oob <- model_data |>
  select(patient_id,
         long_response,
         pfs_time_months,
         disease_progression_to_osimertinib) |>
  mutate(oob_pred = spatial_rf$predicted.oob[, "Yes"])

# Clear environment ----

rm(list = setdiff(ls(), c("spatial_rf",
                          "spatial_rf_oob",
                          "seed",
                          "nobs",
                          "max_ntree",
                          "topvars",
                          "rf_tuning",
                          "rf_tuning_ntree")))

# List with results ----
spatial_rf_list <- list(spatial_rf = spatial_rf,
                        spatial_rf_oob = spatial_rf_oob)

saveRDS(spatial_rf_list,
        file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/6_OSIRESP_long_term_response_signature_spatial.rds")

# PSEUDO-BULK PREDICTION MODEL (SAME GENES) ----

# Load pseudo-bulk data ----

pseudobulk_data <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_pseudobulk.rds")

# Filter out patients with non-evaluable long-term osimertinib response
pseudobulk_data <- pseudobulk_data |>
  filter(long_response != "Non-evaluable") |>
  mutate(long_response = long_response |> factor())

# Data preparation ----

# Special characters in gene names
colnames(pseudobulk_data) <- janitor::make_clean_names(colnames(pseudobulk_data),
                                                       case = "screaming_snake")

# Recover lowercase
pseudobulk_data <- pseudobulk_data |>
  dplyr::rename(patient_id = PATIENT_ID, 
                long_response = LONG_RESPONSE)

# Global parameters ----

topvars <- str_remove_all(topvars, "_stroma|_tumor") |>
  unique()

n_varpro_vars <- topvars |> length()

pseudobulk_data_subset <- pseudobulk_data |>
  select(c("long_response", topvars))

# Number of features
# variables to keep in pre-filtering
n_top <- n_varpro_vars
n_top_sqrt <- n_top |> sqrt() 

nodesize_grid <- floor(c(0.02, 0.05, 0.10, 0.15)*nobs)

mtry_grid <- floor(c(n_top_sqrt/2, n_top_sqrt, n_top_sqrt*1.5, 
                     n_top_sqrt*2, n_top/3, n_top/2, 
                     n_top)) |> unique() |> sort() 

# Split rule
splitrule <- "gini"

# Hyperparameter tuning ----

# 1. mtry and nodesize ----

# Set seed
# set.seed(seed)
rf_tuning_results_pseudo <- rf_tuning(data = pseudobulk_data_subset,
                                      mtry_grid = mtry_grid,
                                      nodesize_grid = nodesize_grid,
                                      max_ntree = max_ntree,
                                      seed = seed,
                                      splitrule = splitrule)

rf_tuning_results_pseudo

# 2. ntree ----

# Set seed
# set.seed(seed)
# ntree
# (using optimal values for mtry and nodesize)
rf_tuning_ntree_results_pseudo <- rf_tuning_ntree(data = pseudobulk_data_subset,
                                                  mtry = rf_tuning_results_pseudo$mtry,
                                                  nodesize = rf_tuning_results_pseudo$nodesize,
                                                  max_ntree = max_ntree,
                                                  seed = seed,
                                                  splitrule = splitrule)

rf_tuning_ntree_results_pseudo

plot(rf_tuning_ntree_results_pseudo$rf_object)

ntree_opt_pseudo <- rf_tuning_ntree_results_pseudo$ntree_opt[, "ntree"] |>
  as.numeric()

# Final RF pseudobulk model ----

# Set seed
# set.seed(seed)
pseudobulk_rf <- imbalanced(long_response ~ ., 
                            data = pseudobulk_data_subset,
                            method = "rfq",
                            splitrule = splitrule,
                            # perf.type = "g.mean",
                            ntree = ntree_opt_pseudo, 
                            mtry = rf_tuning_results_pseudo$mtry,
                            nodesize = rf_tuning_results_pseudo$nodesize,
                            nsplit = 10,
                            # importance = "permute",
                            seed = seed)

saveRDS(pseudobulk_rf,
        file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/6_OSIRESP_long_term_response_signature_pseudobulk.rds")


