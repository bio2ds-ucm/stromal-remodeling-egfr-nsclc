# DESCRIPTION ----
# Required functions for RF gene signature development (script:
# 6_OSIRESP_long_term_response_signature_development.R)

# LIBRARIES ----
library(dplyr)
library(randomForestSRC)

# FUNCTIONS ----

# RF tuning: mtry and nodesize ----
rf_tuning <- function(data,
                      mtry_grid,
                      nodesize_grid,
                      max_ntree,
                      seed,
                      splitrule,
                      case.wt = NULL) {
  tuning_results <- lapply(
    mtry_grid,
    FUN = function(mtry) {
      nodesize_results <- sapply(
        nodesize_grid,
        FUN = function(nodesize) {
          # Random Forest
          rf <- imbalanced(
            long_response ~ .,
            data = data,
            method = "rfq",
            splitrule = splitrule,
            ntree = max_ntree,
            mtry = mtry,
            nodesize = nodesize,
            importance = "none",
            seed = seed,
            case.wt = case.wt
          )
          
          # Out-of-bag AUC
          auc_oob <- get.auc(y = data$long_response,
                             prob = rf$predicted.oob)
          
          # Results
          results <- c("mtry" = mtry,
                       "nodesize" = nodesize,
                       "auc_oob" = auc_oob)
        }
      ) |>
        t()
      
    }
  )
  
  tuning_results <- do.call(rbind, tuning_results)
  
  # Optimal hyperparameter values
  tuning_opt <- tuning_results |>
    as_tibble() |>
    arrange(desc(auc_oob)) |>
    dplyr::slice(1)
  
  tuning_opt
}

# RF tuning: ntree ----
rf_tuning_ntree <- function(data,
                            mtry,
                            nodesize,
                            max_ntree,
                            seed,
                            splitrule,
                            case.wt = NULL) {
  rf_ntree <- imbalanced(
    long_response ~ .,
    data = data,
    method = "rfq",
    splitrule = splitrule,
    ntree = max_ntree,
    mtry = mtry,
    nodesize = nodesize,
    nsplit = 10,
    block.size = 10,
    importance = "none",
    seed = seed,
    case.wt = case.wt
  )
  
  ntree_grid <- seq(10, max_ntree, 10)
  
  auc_oob <- sapply(
    ntree_grid,
    FUN = function(ntree) {
      predict_obj <- predict(rf_ntree, get.tree = 1:ntree)
      
      auc_oob <- get.auc(y = data$long_response,
                         prob = predict_obj$predicted.oob)
      
      results <- c("ntree" = ntree, "auc_oob" = auc_oob)
    }
  ) |>
    t() |>
    as_tibble()
  
  # 1SE rule for determining ntree value
  #     maximum auc value (auc_max)
  auc_max <- auc_oob$auc_oob |> max()
  #     standard deviation of auc values (auc_sd)
  auc_sd <- auc_oob$auc_oob |> sd()
  #     auc_max - auc_sd
  auc_1sd <- auc_max - auc_sd
  
  auc_oob <- auc_oob |>
    mutate(in_1sd_range = (auc_oob >= auc_1sd))
  
  # Minimum ntree with auc >= (auc_max - auc_sd)
  ntree_opt <- auc_oob |>
    arrange(desc(in_1sd_range), ntree) |>
    dplyr::slice(1)
  
  # Return
  list(
    rf_object = rf_ntree,
    ntree_opt = ntree_opt,
    auc_max = auc_max,
    auc_1sd = auc_1sd
  )
}
