# DESCRIPTION ----
# Required functions for ENET logistic regression analyses (scripts:
# 3_OSIRESP_long_term_response_biomarkers_tumor.R,
# 4_OSIRESP_long_term_biomarkers_stroma.R, and
# 5_OSIRESP_long_term_response_biomarkers_joint.R)

# LIBRARIES ----
library(dplyr)
library(furrr)
library(SIS)
library(glmnet)
library(caret)
library(pROC)

# foldid_cv
# creates stratified folds for cross-validation
foldid_cv <- function(y) {
  foldid_list <- createFolds(y |> factor(), k = nfolds)
  foldid_vector <- lapply(
    1:length(foldid_list),
    FUN = function(x) {
      index <- foldid_list[[x]]
      foldid <- rep(x, length(index))
      
      tibble(index = index, foldid = foldid)
    }
  ) |>
    bind_rows() |>
    arrange(index) |>
    select(foldid) |>
    pull()
  
  list(foldid_list = foldid_list, foldid_vector = foldid_vector)
}


# cv_glmnet_alpha

# applies cv.glmnet for each value in alpha grid
# uses family = 'binomial' and type.measure = 'auc'
# returns alpha and lambda that maximize roc auc
cv_glmnet_alpha <- function(x,
                            y,
                            alpha_grid,
                            nlambda,
                            penalty.factor = rep(1, ncol(x)),
                            foldid,
                            nfolds) {
  # cv.glmnet for every value of alpha
  cv_alpha <- lapply(
    alpha_grid,
    FUN = function(alpha) {
      cv <- cv.glmnet(
        x = x,
        y = y,
        alpha = alpha,
        nlambda = nlambda,
        family = 'binomial',
        type.measure = 'auc',
        penalty.factor = penalty.factor,
        nfolds = nfolds,
        foldid = foldid
      )
      cv_results <- data.frame(alpha = alpha,
                               lambda = cv$lambda,
                               auc = cv$cvm)
      cv_results
    }
  ) |>
    bind_rows()
  
  cv_alpha
}

# SIS + elastic net logistic regression model
sis_enet <- function(x_matrix,
                     y_vector,
                     unpen_vars = NULL,
                     unpen_penalties = NULL) {
  # Sure Independence Screening ----
  d_sis <- nrow(x_matrix) - 1
  
  sis <- SIS(
    x_matrix,
    y_vector,
    family = "binomial",
    standardize = TRUE,
    iter = FALSE,
    nsis = d_sis
  )
  
  sis_ind <- sis$sis.ix0
  
  # Penalty factor ----
  
  # Inf for excluded variables
  # 1 for included and penalized (default penalty factor)
  # 0 for non-penalized variables
  
  # Penalty factor initialization
  penalty_factor <- rep(1, ncol(x_matrix))
  
  # Genes filtered out through SIS method
  penalty_factor[setdiff(1:ncol(x_matrix), sis_ind)] <- Inf
  
  # Unpenalized/less penalized predictors
  if (!is.null(unpen_vars)) {
    for (ind in 1:length(unpen_vars)) {
      penalty_factor[which(colnames(x_matrix) == unpen_vars[ind])] <- unpen_penalties[ind]
    }
  }
  
  # Hyperparameter tuning
  
  # k-fold cross-validation for hyperparameter tuning
  # criterion: maximization of ROC AUC
  cv <- lapply(
    1:cvrep_tuning,
    FUN = function(rep) {
      foldid <- foldid_cv(y = y_vector)
      
      cv_glmnet_alpha(
        x = x_matrix,
        y = y_vector,
        alpha_grid = alpha_grid,
        nlambda = nlambda,
        foldid = foldid$foldid_vector,
        nfolds = nfolds,
        penalty.factor = penalty_factor
      )
    }
  ) |>
    bind_rows()
  
  # For each alpha, get q3 from lambda.1min distribution
  lambda_perc <- cv |>
    group_by(alpha) |>
    summarise_at(vars(lambda), list(
      lambda_perc = function(x) {
        quantile(x, probs = 0.75)
      }
    ))
  
  # Join alpha, lambda and auc values
  cv_opt <- cv |>
    left_join(lambda_perc, by = join_by(alpha)) |>
    mutate(lambda_diff = abs(lambda - lambda_perc)) |>
    group_by(alpha) |>
    arrange(lambda_diff)  |>
    slice(1) |>
    ungroup() |>
    arrange(desc(auc)) |>
    select(alpha, lambda) |>
    slice(1)
  
  # Optimum values for alpha and lambda
  alpha_opt <- cv_opt |>
    select(alpha) |>
    pull()
  
  lambda_opt <- cv_opt |>
    select(lambda) |>
    pull()
  
  # Model
  model <- glmnet(
    x = x_matrix,
    y = y_vector,
    alpha = alpha_opt,
    nlambda = nlambda,
    family = "binomial",
    penalty.factor = penalty_factor
  )
  
  # Model coefficients
  model_coef <- predict(model, type = "coefficients", s = lambda_opt) |>
    as.matrix() |>
    as.data.frame()
  
  # Results
  list(
    model = model,
    alpha_opt = alpha_opt,
    lambda_opt = lambda_opt,
    model_coef = model_coef
  )
}

# Predicted probabilities
predictions <- function(model, lambda_opt, newx) {
  # Predictions
  predict(model,
          type = "response",
          s = lambda_opt,
          newx = newx) |>
    as.vector()
  
}

# ROC AUC estimates
auc_estimates <- function(model_data, y_vector, pred) {
  # Recode ECOG 0-1 in case of being numeric
  if (is.numeric(model_data$ecog_ps_dummie)) {
    model_data <- model_data |>
      mutate(ecog_ps_dummie = case_when(ecog_ps_dummie == 0 ~ "< 2", .default = "\u2265 2"))
  }
  
  # Global AUC
  auc <- roc(response = y_vector,
             predictor = pred,
             direction = "<")$auc
  
  # AUC typical EGFR mutations
  typical_ind <- model_data$egfr_mutation_type != "Other"
  
  if ((y_vector[typical_ind] |> table() |> length()) > 1) {
    auc_typical <- roc(response = y_vector[typical_ind],
                       predictor = pred[typical_ind],
                       direction = "<")$auc
  }
  else{
    auc_typical <- NA
  }
  
  # AUC ECOG 0-1
  ecog_ind <- model_data$ecog_ps_dummie == "< 2"
  
  if ((y_vector[ecog_ind] |> table() |> length()) > 1) {
    auc_ecog <- roc(response = y_vector[ecog_ind],
                    predictor = pred[ecog_ind],
                    direction = "<")$auc
  }
  else{
    auc_ecog <- NA
  }
  
  # AUC ECOG 0-1 typical EGFR mutations
  typical_ecog_ind <- model_data$egfr_mutation_type != "Other" &
    model_data$ecog_ps_dummie == "< 2"
  
  if ((y_vector[typical_ecog_ind] |> table() |> length()) > 1) {
    auc_typical_ecog <- roc(response = y_vector[typical_ecog_ind],
                            predictor = pred[typical_ecog_ind],
                            direction = "<")$auc
  }
  else{
    auc_typical_ecog <- NA
  }
  
  # Results
  list(
    auc = auc,
    auc_typical = auc_typical,
    auc_ecog = auc_ecog,
    auc_typical_ecog = auc_typical_ecog
  )
}

# Internal validation (repeated stratified k-fold cross-validation)
repcv <- function(model_data,
                  x_matrix,
                  y_vector,
                  unpen_vars = NULL,
                  unpen_penalties = NULL) {
  future::plan(multisession, workers = ncores)
  
  # Iterate over repetitions of 5-fold cross-validation
  result <- future_map(1:cvrep, function(rep) {
    # Stratified folds
    folds <- createFolds(y = y_vector |> as.factor(), k = nfolds)
    
    # Iterate over test folds
    kfold_cv <- lapply(
      1:nfolds,
      FUN = function(fold) {
        # Observations selected in test fold
        testfold_ind <- folds[[fold]]
        
        # Subset train data
        x_matrix_trainfolds <- x_matrix[-testfold_ind, ]
        y_vector_trainfolds <- y_vector[-testfold_ind]
        
        # Subset test data
        x_matrix_testfold <- x_matrix[testfold_ind, ]
        y_vector_testfold <- y_vector[testfold_ind]
        
        model_results_fold <- sis_enet(
          x_matrix = x_matrix_trainfolds,
          y_vector = y_vector_trainfolds,
          unpen_vars = unpen_vars,
          unpen_penalties = unpen_penalties
        )
        
        # Predictions in train folds
        pred_trainfolds <- predictions(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          newx = x_matrix_trainfolds
        )
        
        # Prediccions in test fold
        pred_testfold <- predictions(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          newx = x_matrix_testfold
        )
        
        # ROC AUC
        auc_trainfolds <- auc_estimates(model_data = model_data[-testfold_ind, ],
                                        y_vector = y_vector_trainfolds,
                                        pred = pred_trainfolds)
        
        auc_testfold <- auc_estimates(model_data = model_data[testfold_ind, ],
                                      y_vector = y_vector_testfold,
                                      pred = pred_testfold)
        
        # Results (predictions)
        nrow_results <- model_data |> nrow()
        pred_df <- data.frame(
          rep = rep(rep, nrow_results),
          test_fold = rep(fold, nrow_results),
          auc_trainfolds = rep(auc_trainfolds$auc, nrow_results),
          auc_testfold = rep(auc_testfold$auc, nrow_results),
          auc_typical_trainfolds = rep(auc_trainfolds$auc_typical, nrow_results),
          auc_typical_testfold = rep(auc_testfold$auc_typical, nrow_results),
          auc_ecog_trainfolds = rep(auc_trainfolds$auc_ecog, nrow_results),
          auc_ecog_testfold = rep(auc_testfold$auc_ecog, nrow_results),
          auc_ecog_typical_trainfolds = rep(auc_trainfolds$auc_typical_ecog, nrow_results),
          auc_ecog_typical_testfold = rep(auc_testfold$auc_typical_ecog, nrow_results),
          patient_id = model_data$patient_id,
          pred_trainfolds = NA,
          pred_testfold = NA
        )
        
        pred_df[setdiff(1:nrow_results, testfold_ind), "pred_trainfolds"] <- pred_trainfolds
        pred_df[testfold_ind, "pred_testfold"] <- pred_testfold
        
        # Results (coefficients)
        n_nonull_coef <- model_results_fold$model_coef |>
          as.data.frame() |>
          filter(s1 != 0) |>
          nrow()
        
        # Intercept must not be considered
        n_nonull_coef <- n_nonull_coef - 1
        
        model_coef_df <- model_results_fold$model_coef |>
          t() |>
          as.data.frame() |>
          mutate(
            rep = rep,
            testfold = fold,
            alpha = model_results_fold$alpha_opt,
            lambda = model_results_fold$lambda_opt,
            n_nonull_coef = n_nonull_coef
          )
        
        # Results
        list(prediccions = pred_df, coefficients = model_coef_df)
      }
    )
    
  }, .options = furrr_options(seed = TRUE))
  
  plan(sequential)
  
  result
}