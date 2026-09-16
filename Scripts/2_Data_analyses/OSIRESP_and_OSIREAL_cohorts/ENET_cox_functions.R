# DESCRIPTION ----
# Required functions for ENET Cox regression analyses (scripts:
# 1_OSIRESP_resistance_biomarkers_tumor.R and
# 2_OSIRESP_resistance_biomarkers_stroma.R)

# LIBRARIES ----
library(dplyr)
library(furrr)
library(SIS)
library(glmnet)
library(caret)
library(timeROC)

# FUNCTIONS ----

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

# cv_glmnet_alpha_cox

# applies cv.glmnet for each value in alpha grid
# uses family = "cox" and type.measure = "C"
# returns alpha and lambda that maximize C index
cv_glmnet_alpha_cox <- function(x,
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
        family = "cox",
        type.measure = "C",
        penalty.factor = penalty.factor,
        nfolds = nfolds,
        foldid = foldid
      )
      
      cv_results <- data.frame(alpha = alpha,
                               lambda = cv$lambda,
                               c_index = cv$cvm)
      cv_results
    }
  ) |>
    bind_rows()
  
  cv_alpha
}

# SIS + ENET Cox regression
sis_enet_cox <- function(x_matrix,
                         y_surv,
                         unpen_vars = NULL,
                         unpen_penalties = NULL) {
  # SIS ----
  d_sis <- nrow(x_matrix) - 1
  
  # For family = "cox", y should be an object of class Surv, as provided by the
  # function Surv() in the package survival.
  sis <- SIS(
    x_matrix,
    y_surv,
    family = "cox",
    penalty = "lasso",
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
  
  # Unpenalized predictors
  if (!is.null(unpen_vars)) {
    for (ind in 1:length(unpen_vars)) {
      penalty_factor[which(colnames(x_matrix) == unpen_vars[ind])] <- unpen_penalties[ind]
    }
  }
  
  # Hyperparameter tuning
  
  # k-fold cross-validation for hyperparameter tuning
  # criterion: maximization of C index
  cv <- lapply(
    1:cvrep_tuning,
    FUN = function(rep) {
      foldid <- foldid_cv(y = y_surv[, 2])
      
      cv_glmnet_alpha_cox(
        x = x_matrix,
        y = y_surv,
        alpha_grid = alpha_grid,
        nlambda = nlambda,
        foldid = foldid$foldid_vector,
        nfolds = nfolds,
        penalty.factor = penalty_factor
      )
    }
  ) |>
    bind_rows()
  
  # For each alpha, get percentile from lambda.1min distribution
  lambda_perc <- cv |>
    group_by(alpha) |>
    summarise_at(vars(lambda), list(
      lambda_perc = function(x) {
        quantile(x, probs = 0.75)
      }
    ))
  
  # Join alpha, lambda and C index values
  cv_opt <- cv |>
    left_join(lambda_perc, by = join_by(alpha)) |>
    mutate(lambda_diff = abs(lambda - lambda_perc)) |>
    group_by(alpha) |>
    arrange(lambda_diff)  |>
    slice(1) |>
    ungroup() |>
    arrange(desc(c_index)) |>
    select(alpha, lambda) |>
    slice(1)
  
  # Optimal values for alpha and lambda
  alpha_opt <- cv_opt |>
    select(alpha) |>
    pull()
  
  lambda_opt <- cv_opt |>
    select(lambda) |>
    pull()
  
  # Model
  model <- glmnet(
    x = x_matrix,
    y = y_surv,
    alpha = alpha_opt,
    nlambda = nlambda,
    family = "cox",
    penalty.factor = penalty_factor
  )
  
  # Model coefficients
  model_coef <- predict(model, type = "coefficients", s = lambda_opt) |>
    as.matrix() |>
    as.data.frame()
  colnames(model_coef) <- "s1"
  
  # Results
  list(
    model = model,
    alpha_opt = alpha_opt,
    lambda_opt = lambda_opt,
    model_coef = model_coef
  )
}

# Linear predictor
predictions_cox <- function(model, lambda_opt, newx) {
  # Predictions
  predict(model,
          type = "link",
          s = lambda_opt,
          newx = newx) |>
    as.vector()
}

# C index estimates
c_index_estimates <- function(model,
                              lambda_opt,
                              model_data,
                              x_matrix,
                              y_surv) {
  # Global C index
  c_index <- assess.glmnet(
    model,
    newx = x_matrix,
    newy = y_surv,
    family = "cox",
    s = lambda_opt
  )$C
  
  
  # C index typical EGFR mutations
  typical_ind <- model_data$egfr_mutation_type != "Other"
  
  c_index_typical <- assess.glmnet(
    model,
    newx = x_matrix[typical_ind, ],
    newy = y_surv[typical_ind, ],
    family = "cox",
    s = lambda_opt
  )$C
  
  # C index ECOG 0-1
  ecog_ind <- model_data$ecog_ps_dummie == "< 2"
  
  c_index_ecog <- assess.glmnet(
    model,
    newx = x_matrix[ecog_ind, ],
    newy = y_surv[ecog_ind, ],
    family = "cox",
    s = lambda_opt
  )$C
  
  # C index ECOG 0-1 typical EGFR mutations
  typical_ecog_ind <- model_data$egfr_mutation_type != "Other" &
    model_data$ecog_ps_dummie == "< 2"
  
  c_index_typical_ecog <- assess.glmnet(
    model,
    newx = x_matrix[typical_ecog_ind, ],
    newy = y_surv[typical_ecog_ind, ],
    family = "cox",
    s = lambda_opt
  )$C
  
  # Results
  list(
    c_index = c_index,
    c_index_typical = c_index_typical,
    c_index_ecog = c_index_ecog,
    c_index_typical_ecog = c_index_typical_ecog
  )
}


# Time dependent ROC AUC estimates
timeroc_estimates <- function(model_data,
                              y_times,
                              y_events,
                              marker,
                              cause,
                              times,
                              iid = FALSE) {
  # Global
  timeroc <- timeROC(
    T = y_times,
    delta = y_events,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = times,
    ROC = TRUE,
    iid = iid
  )
  
  
  # Typical EGFR mutations
  typical_ind <- model_data$egfr_mutation_type != "Other"
  
  timeroc_typical <- timeROC(
    T = y_times[typical_ind],
    delta = y_events[typical_ind],
    marker = marker[typical_ind],
    cause = 1,
    weighting = "marginal",
    times = times,
    ROC = TRUE,
    iid = iid
  )
  
  # ECOG 0-1
  ecog_ind <- model_data$ecog_ps_dummie == "< 2"
  
  timeroc_ecog <- timeROC(
    T = y_times[ecog_ind],
    delta = y_events[ecog_ind],
    marker = marker[ecog_ind],
    cause = 1,
    weighting = "marginal",
    times = times,
    ROC = TRUE,
    iid = iid
  )
  
  # ECOG 0-1 typical EGFR mutations
  typical_ecog_ind <- model_data$egfr_mutation_type != "Other" &
    model_data$ecog_ps_dummie == "< 2"
  
  timeroc_typical_ecog <- timeROC(
    T = y_times[typical_ecog_ind],
    delta = y_events[typical_ecog_ind],
    marker = marker[typical_ecog_ind],
    cause = 1,
    weighting = "marginal",
    times = times,
    ROC = TRUE,
    iid = iid
  )
  
  # Results
  list(
    timeroc = timeroc,
    timeroc_typical = timeroc_typical,
    timeroc_ecog = timeroc_ecog,
    timeroc_typical_ecog = timeroc_typical_ecog
  )
}

# Predicted survival probabilities at specified timepoints

survfit_prob <- function(model,
                         lambda_opt,
                         x_matrix,
                         y_surv,
                         newx,
                         times) {
  survf <- survfit(
    model,
    s = lambda_opt,
    x = x_matrix,
    y = y_surv,
    newx = newx
  )
  
  survf_times <- summary(survf, times = times)$surv
}

# Internal validation (repeated stratified k-fold cross-validation)
repcv_cox <- function(model_data,
                      x_matrix,
                      y_surv,
                      unpen_vars = NULL,
                      unpen_penalties = NULL) {
  future::plan(multisession, workers = ncores)
  
  # Iterate over repetitions of 5-fold cross-validation
  result <- future_map(1:cvrep, function(rep) {
    # Stratified folds
    folds <- createFolds(y = y_surv[, 2] |> as.factor(), k = nfolds)
    
    # Iterate over test folds
    kfold_cv <- lapply(
      1:nfolds,
      FUN = function(fold) {
        # Observations selected in test fold
        testfold_ind <- folds[[fold]]
        
        # Subset train data
        x_matrix_trainfolds <- x_matrix[-testfold_ind, ]
        y_surv_trainfolds <- y_surv[-testfold_ind]
        x_matrix_trainfolds |> dim()
        y_surv_trainfolds |> dim()
        
        # Subset test data
        x_matrix_testfold <- x_matrix[testfold_ind, ]
        y_surv_testfold <- y_surv[testfold_ind]
        x_matrix_testfold |> dim()
        y_surv_testfold |> dim()
        
        model_results_fold <- sis_enet_cox(
          x_matrix = x_matrix_trainfolds,
          y_surv = y_surv_trainfolds,
          unpen_vars = unpen_vars,
          unpen_penalties = unpen_penalties
        )
        
        # Linear predictor in train folds
        pred_trainfolds <- predictions_cox(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          newx = x_matrix_trainfolds
        )
        
        # Linear predictor in test folds
        pred_testfold <- predictions_cox(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          newx = x_matrix_testfold
        )
        
        # C index
        c_index_trainfolds <- c_index_estimates(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          model_data = model_data[-testfold_ind, ],
          x_matrix = x_matrix_trainfolds,
          y_surv = y_surv_trainfolds
        )
        
        c_index_testfold <- c_index_estimates(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          model_data = model_data[testfold_ind, ],
          x_matrix = x_matrix_testfold,
          y_surv = y_surv_testfold
        )
        
        # Time dependent ROC AUC
        timeroc_trainfolds <- timeroc_estimates(
          model_data = model_data[-testfold_ind, ],
          y_times = y_surv_trainfolds[, 1],
          y_events = y_surv_trainfolds[, 2],
          marker = pred_trainfolds,
          cause = 1,
          times = times
        )
        
        timeroc_testfold <- timeroc_estimates(
          model_data = model_data[testfold_ind, ],
          y_times = y_surv_testfold[, 1],
          y_events = y_surv_testfold[, 2],
          marker = pred_testfold,
          cause = 1,
          times = times
        )
        
        # Predicted survival probabilities at specified timepoints
        surv_prob_trainfolds <- survfit_prob(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          x_matrix = x_matrix_trainfolds,
          y_surv = y_surv_trainfolds,
          newx = x_matrix_trainfolds,
          times = times
        )
        
        surv_prob_testfold <- survfit_prob(
          model = model_results_fold$model,
          lambda_opt = model_results_fold$lambda_opt,
          x_matrix = x_matrix_trainfolds,
          y_surv = y_surv_trainfolds,
          newx = x_matrix_testfold,
          times = times
        )
        
        # Results (predictions)
        nrow_results <- model_data |> nrow()
        pred_df <- data.frame(
          rep = rep(rep, nrow_results),
          test_fold = rep(fold, nrow_results),
          c_index_trainfolds = rep(c_index_trainfolds$c_index, nrow_results),
          c_index_testfold = rep(c_index_testfold$c_index, nrow_results),
          c_index_typical_trainfolds = rep(c_index_trainfolds$c_index_typical, nrow_results),
          c_index_typical_testfold = rep(c_index_testfold$c_index_typical, nrow_results),
          c_index_ecog_trainfolds = rep(c_index_trainfolds$c_index_ecog, nrow_results),
          c_index_ecog_testfold = rep(c_index_testfold$c_index_ecog, nrow_results),
          c_index_ecog_typical_trainfolds = rep(c_index_trainfolds$c_index_typical_ecog, nrow_results),
          c_index_ecog_typical_testfold = rep(c_index_testfold$c_index_typical_ecog, nrow_results),
          timeroc10_trainfolds = rep(timeroc_trainfolds$timeroc$AUC[1], nrow_results),
          timeroc10_testfold = rep(timeroc_testfold$timeroc$AUC[1], nrow_results),
          timeroc10_typical_trainfolds = rep(timeroc_trainfolds$timeroc_typical$AUC[1], nrow_results),
          timeroc10_typical_testfold = rep(timeroc_testfold$timeroc_typical$AUC[1], nrow_results),
          timeroc10_ecog_trainfolds = rep(timeroc_trainfolds$timeroc_ecog$AUC[1], nrow_results),
          timeroc10_ecog_testfold = rep(timeroc_testfold$timeroc_ecog$AUC[1], nrow_results),
          timeroc10_ecog_typical_trainfolds = rep(
            timeroc_trainfolds$timeroc_typical_ecog$AUC[1],
            nrow_results
          ),
          timeroc10_ecog_typical_testfold = rep(
            timeroc_testfold$timeroc_typical_ecog$AUC[1],
            nrow_results
          ),
          timeroc24_trainfolds = rep(timeroc_trainfolds$timeroc$AUC[2], nrow_results),
          timeroc24_testfold = rep(timeroc_testfold$timeroc$AUC[2], nrow_results),
          timeroc24_typical_trainfolds = rep(timeroc_trainfolds$timeroc_typical$AUC[2], nrow_results),
          timeroc24_typical_testfold = rep(timeroc_testfold$timeroc_typical$AUC[2], nrow_results),
          timeroc24_ecog_trainfolds = rep(timeroc_trainfolds$timeroc_ecog$AUC[2], nrow_results),
          timeroc24_ecog_testfold = rep(timeroc_testfold$timeroc_ecog$AUC[2], nrow_results),
          timeroc24_ecog_typical_trainfolds = rep(
            timeroc_trainfolds$timeroc_typical_ecog$AUC[2],
            nrow_results
          ),
          timeroc24_ecog_typical_testfold = rep(
            timeroc_testfold$timeroc_typical_ecog$AUC[2],
            nrow_results
          ),
          patient_id = model_data$patient_id,
          pred_trainfolds = NA,
          pred_testfold = NA,
          surv_prob10_trainfolds = NA,
          surv_prob24_trainfolds = NA,
          surv_prob10_testfold = NA,
          surv_prob24_testfold = NA
        )
        
        pred_df[setdiff(1:nrow_results, testfold_ind), "pred_trainfolds"] <- pred_trainfolds
        pred_df[testfold_ind, "pred_testfold"] <- pred_testfold
        
        pred_df[setdiff(1:nrow_results, testfold_ind), "surv_prob10_trainfolds"] <- surv_prob_trainfolds[1, ]
        pred_df[setdiff(1:nrow_results, testfold_ind), "surv_prob24_trainfolds"] <- surv_prob_trainfolds[2, ]
        
        pred_df[testfold_ind, "surv_prob10_testfold"] <- surv_prob_testfold[1, ]
        pred_df[testfold_ind, "surv_prob24_testfold"] <- surv_prob_testfold[2, ]
        
        
        # Results (coefficients)
        n_nonull_coef <- model_results_fold$model_coef |>
          as.data.frame() |>
          filter(s1 != 0) |>
          nrow()
        
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
