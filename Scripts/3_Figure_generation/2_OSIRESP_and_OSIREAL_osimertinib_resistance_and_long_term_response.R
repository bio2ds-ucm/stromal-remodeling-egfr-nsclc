# DESCRIPTION ----
# Figures 2, S1, and S2


# LIBRARIES ----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(ggridges)
library(tidytext)
library(randomForestSRC)
library(pROC)
library(survival)
library(survminer)
library(png)
library(cowplot)
library(gridExtra)
library(ggplotify)

# COMPARTMENT-SPECIFIC OSIMERTINIB RESISTANCE BIOMARKERS -----

# 1. TUMOR COMPARTMENT ----

# Load results
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/1_OSIRESP_resistance_biomarkers_tumor.RData"
)

# Save object with results using a different name
model_tumor <- model_results

# Prediction results
repcv_pred_tumor <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$prediccions
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor")

# Coefficients per gene
repcv_coef_tumor <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$coefficients
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor")

# Check dimensions. Total number of elastic net models is 100 rep x 5 fold = 500
repcv_coef_tumor |> dim()

# Results per gene:
# - number of models where selected
# - proportion of the total of 500 models
# - minimum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - maximum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - mean coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - median coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
repcv_coef_summary_tumor <- repcv_coef_tumor |>
  select(-rep, -testfold, -alpha, -lambda, -n_nonull_coef, -compartment) |>
  sapply(
    FUN = function(x) {
      nonzero <- sum(x != 0)
      nonzero_prop <- nonzero / length(x)
      x_nonzero <- x[x != 0]
      min <- min(x_nonzero)
      mean <- mean(x_nonzero)
      median <- median(x_nonzero)
      max <- max(x_nonzero)
      results <- c(
        "gene" = names(x),
        "nonzero_values" = nonzero,
        "nonzero_prop" = nonzero_prop,
        "min_coef" = min,
        "max_coef" = max,
        "mean_coef" = mean,
        "median_coef" = median
      )
      results
    }
  ) |>
  t()

# Gene names
gene_names <- repcv_coef_summary_tumor |>
  rownames()

# Replace with 0 values
# Add gene names
# Add compartmnent
# Add overall sign of association with the outcome
# Set order of columns
repcv_coef_summary_tumor <- repcv_coef_summary_tumor |>
  as_tibble() |>
  mutate_all(
    .funs = function(x) {
      replace(x, !is.finite(x), 0)
    }
  ) |>
  mutate(
    gene = gene_names,
    compartment = "tumor",
    overall_sign = case_when(
      max_coef <= 0 ~ "negative association with osimertinib resistance",
      min_coef >= 0 ~ "positive association with osimertinib resistance",
      .default = "non-conclusive"
    )
  ) |>
  select(
    gene,
    compartment,
    nonzero_values,
    nonzero_prop,
    min_coef,
    max_coef,
    mean_coef,
    median_coef,
    overall_sign
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "model_tumor",
    "repcv_pred_tumor",
    "repcv_coef_tumor",
    "repcv_coef_summary_tumor"
  )
))

# 2. STROMAL COMPARTMENT ----

# Load results
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/2_OSIRESP_resistance_biomarkers_stroma.RData"
)

# Save object with results using a different name
model_stroma <- model_results

# Prediction results
repcv_pred_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$prediccions
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "stroma")

# Coefficients per gene
repcv_coef_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$coefficients
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "stroma")

# Check dimensions. Total number of elastic net models is 100 rep x 5 fold = 500
repcv_coef_stroma |> dim()

# Results per gene:
# - number of models where selected
# - proportion of the total of 500 models
# - minimum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - maximum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - mean coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - median coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
repcv_coef_summary_stroma <- repcv_coef_stroma |>
  select(-rep, -testfold, -alpha, -lambda, -n_nonull_coef, -compartment) |>
  sapply(
    FUN = function(x) {
      nonzero <- sum(x != 0)
      nonzero_prop <- nonzero / length(x)
      x_nonzero <- x[x != 0]
      min <- min(x_nonzero)
      mean <- mean(x_nonzero)
      median <- median(x_nonzero)
      max <- max(x_nonzero)
      results <- c(
        "gene" = names(x),
        "nonzero_values" = nonzero,
        "nonzero_prop" = nonzero_prop,
        "min_coef" = min,
        "max_coef" = max,
        "mean_coef" = mean,
        "median_coef" = median
      )
      results
    }
  ) |>
  t()

# Gene names
gene_names <- repcv_coef_summary_stroma |>
  rownames()

# Replace with 0 values
# Add gene names
# Add compartmnent
# Add overall sign of association with the outcome
# Set order of columns
repcv_coef_summary_stroma <- repcv_coef_summary_stroma |>
  as_tibble() |>
  mutate_all(
    .funs = function(x) {
      replace(x, !is.finite(x), 0)
    }
  ) |>
  mutate(
    gene = gene_names,
    compartment = "stroma",
    overall_sign = case_when(
      max_coef <= 0 ~ "negative association with osimertinib resistance",
      min_coef >= 0 ~ "positive association with osimertinib resistance",
      .default = "non-conclusive"
    )
  ) |>
  select(
    gene,
    compartment,
    nonzero_values,
    nonzero_prop,
    min_coef,
    max_coef,
    mean_coef,
    median_coef,
    overall_sign
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "model_tumor",
    "repcv_pred_tumor",
    "repcv_coef_tumor",
    "repcv_coef_summary_tumor",
    "model_stroma",
    "repcv_pred_stroma",
    "repcv_coef_stroma",
    "repcv_coef_summary_stroma"
  )
))

# 3. BIND RESULTS ----

# For graphical representation
# - prediction results
repcv_pred <- list(repcv_pred_tumor, repcv_pred_stroma) |>
  bind_rows() |>
  mutate(compartment = compartment |>
           factor(levels = c("tumor", "stroma"), ordered = TRUE))

# - coefficients
repcv_coef <- list(repcv_coef_tumor, repcv_coef_stroma) |>
  bind_rows() |>
  mutate(compartment = compartment |>
           factor(levels = c("tumor", "stroma"), ordered = TRUE))

# - coefficients summary
repcv_coef_summary <- list(repcv_coef_summary_tumor, repcv_coef_summary_stroma) |>
  bind_rows() |>
  mutate(
    gene = str_remove_all(gene, "`"),
    compartment = compartment |>
      factor(levels = c("tumor", "stroma"), ordered = TRUE),
    overall_sign = overall_sign |>
      factor(
        levels = c(
          "positive association with osimertinib resistance",
          "negative association with osimertinib resistance",
          "non-conclusive"
        ),
        ordered = TRUE
      )
  )

# 4. FIGURES ----

# Lollipot plot: robust biomarkers ----

gene_inclusion_plot <- ggplot(repcv_coef_summary |>
                                filter(nonzero_prop >= 0.40)) +
  geom_segment(
    aes(
      x = reorder_within(gene, nonzero_prop, compartment),
      xend = reorder_within(gene, nonzero_prop, compartment),
      y = nonzero_prop,
      yend = 0,
      color = overall_sign
    ),
    size = 0.05
  ) +
  geom_point(aes(
    x = reorder_within(gene, nonzero_prop, compartment),
    y = nonzero_prop,
    color = overall_sign
  ),
  size = 1) +
  scale_x_reordered() +
  scale_color_manual(
    values = c(
      "negative association with osimertinib resistance" = "darkviolet",
      "positive association with osimertinib resistance" = "brown1"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent
  ) +
  labs(x = "", y = "selection rate", color = "") +
  coord_flip() +
  facet_wrap(~ compartment, scales = "free") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 5),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 6),
    axis.ticks.y = element_line(size = 0.1),
    axis.ticks.x = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    legend.box.margin = margin(-5, 0, -5, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing = unit(0, "pt"),
    strip.background = element_rect(size = 0.25)
  )

# Ridge plot: coefficients/effects across all repetitions ----

# Genes to subset
# (those with a percentage of inclusion >= 40%)
genes_to_subset_tumor <- repcv_coef_summary_tumor |>
  filter(nonzero_prop >= 0.40) |>
  select(gene) |>
  pull()

genes_to_subset_tumor |>
  length()

genes_to_subset_stroma <- repcv_coef_summary_stroma |>
  filter(nonzero_prop >= 0.40) |>
  select(gene) |>
  pull()

genes_to_subset_stroma |>
  length()

# Data in long format
# Filter selected genes based on percentage of inclusion
# Filter null values for regression coefficients
# Remove "`" from gene names (for plot)

# Check colnames (for pivot_longer)
repcv_coef |> colnames() |> head(n = 20)
repcv_coef |> colnames() |> tail(n = 20)

repcv_coef_long <- repcv_coef |>
  pivot_longer(ACTA2:TNFSF12, names_to = "gene", values_to = "coef") |>
  filter((compartment == "tumor" &
            gene %in% genes_to_subset_tumor) |
           (compartment == "stroma" &
              gene %in% genes_to_subset_stroma)
  ) |>
  filter(coef != 0) |>
  mutate(gene = str_remove_all(gene, "`"))

# Add mean coefficient per gene
repcv_coef_long <- repcv_coef_long |>
  left_join(repcv_coef_summary |>
              select(gene, compartment, mean_coef),
            by = join_by(gene, compartment))


# Ridge plot of non-zero regression coefficient values per gen
gene_coef_plot <- ggplot(repcv_coef_long, aes(
  x = coef,
  y = reorder(gene, mean_coef),
  fill = mean_coef
)) +
  geom_vline(xintercept = 0, lwd = 0.05) +
  geom_density_ridges(
    linewidth = 0.05,
    color = "black",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0, height = 0),
    point_shape = "|",
    point_size = 0.5,
    alpha = 0.7
  ) +
  geom_density_ridges(linewidth = 0.05, color = "black") +
  scale_fill_gradientn(
    colours = c("darkviolet", "white", "brown1"),
    values = scales::rescale(c(-1, 0, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  facet_wrap(~ compartment, scales = "free") +
  labs(x = "regression coefficient", y = "", fill = "mean value") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(
    axis.title.x = element_text(hjust = 0.5),
    line = element_line(linewidth = 0.05, color = "black"),
    text = element_text(size = 6),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    strip.background = element_rect(color = "black", linewidth = 0.25),
    legend.key.height = unit(0.01, "npc")
  )

# AUC values ----

repcv_auc <- repcv_pred |>
  select(compartment,
         rep,
         test_fold,
         timeroc10_testfold,
         timeroc24_testfold) |>
  unique()

# Check dimensions
repcv_auc |>
  group_by(compartment) |>
  summarise(n = n())

# Long format
repcv_auc <- repcv_auc |>
  pivot_longer(
    cols = timeroc10_testfold:timeroc24_testfold,
    names_to = "time",
    values_to = "auc_testfold"
  )

alpha <- 0.05

# Get mean values per repetition
repcv_auc_mean <- repcv_auc |>
  group_by(compartment, rep, time) |>
  summarise(auc_testfold = mean(auc_testfold))

# Check dimensions
repcv_auc_mean |>
  group_by(compartment) |>
  summarise(n = n())

# Get point estimate and confidence intervals
repcv_auc_stats <- repcv_auc_mean |>
  group_by(compartment, time) |>
  summarise(
    mean = mean(auc_testfold),
    ci_low = quantile(auc_testfold, probs = alpha / 2),
    ci_upp = quantile(auc_testfold, probs = 1 - alpha / 2),
    label = paste0(
      format(round(mean, digits = 2), nsmall = 2),
      " (",
      format(round(ci_low, digits = 2), nsmall = 2),
      ", ",
      format(round(ci_upp, digits = 2), nsmall = 2),
      ")"
    )
  )

repcv_auc_stats

dodge_width <- 0.1

auc_plot <- ggplot(data = repcv_auc_mean, aes(x = compartment, y = auc_testfold, fill = compartment)) +
  stat_halfeye(point_interval = NULL, alpha = 0.8) +
  geom_jitter(
    data = repcv_auc_mean,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(compartment) - dodge_width,
      y = auc_testfold,
      color = compartment
    ),
    height = 0,
    width = 0.025,
    # alpha = 0.5,
    size = 0.25
  ) +
  geom_point(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(x = compartment, y = mean),
    size = 1
  ) +
  geom_errorbar(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(
      x = compartment,
      y = mean,
      ymin = ci_low,
      ymax = ci_upp
    ),
    width = 0.1,
    linewidth = 0.25
  ) +
  facet_wrap(~ time, labeller = labeller(
    time = c(
      timeroc10_testfold = "resistance at 10 months",
      timeroc24_testfold = "resistance at 24 months"
    ),
  )) +
  scale_y_continuous(limits = c(0.5, 0.9),
                     breaks = seq(0.5, 0.9, by = 0.05)) +
  scale_x_discrete(expand = c(0.25, 0)) +
  # scale_color_manual(values = c("tumor" = "#08CC3C",
  #                               "stroma" = "#121DEE")) +
  # scale_fill_manual(values = c("tumor" = "#08CC3C",
  #                              "stroma" = "#121DEE")) +
  # scale_color_manual(values = c(
  #   "tumor" = "#2EAD73",
  #   "stroma" = "#121DEE"
  # )) +
  # scale_fill_manual(values = c(
  #   "tumor" = "#2EAD73",
  #   "stroma" = "#121DEE"
  # )) +
  scale_color_manual(values = c("tumor" = "#66DEAA", "stroma" = "#4A5FFF")) +
  scale_fill_manual(values = c("tumor" = "#66DEAA", "stroma" = "#4A5FFF")) +
  geom_label(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(x = compartment, y = mean, label = label),
    hjust = -0.15,
    # vjust = 0.2,
    show.legend = FALSE,
    size = 1.5,
    label.size = 0.1
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0),
    text = element_text(size = 6),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25)
  ) +
  labs(x = "", y = "tAUC estimates from repeated cross-validation")

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot"
  )
))

# COMPARTMENT-SPECIFIC LONG-TERM OSIMERTINIB RESPONSE BIOMARKERS ----

# 1. TUMOR COMPARTMENT ----

# Load results
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/3_OSIRESP_long_term_response_biomarkers_tumor.RData"
)

# Save object with results using a different name
model_tumor <- model_results

# Prediction results
repcv_pred_tumor <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$prediccions
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor")

# Coefficients per gene
repcv_coef_tumor <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$coefficients
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor")

# Check dimensions. Total number of elastic net models is 100 rep x 5 fold = 500
repcv_coef_tumor |> dim()

# Results per gene:
# - number of models where selected
# - proportion of the total of 500 models
# - minimum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - maximum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - mean coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - median coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
repcv_coef_summary_tumor <- repcv_coef_tumor |>
  select(-`(Intercept)`,-rep,-testfold,-alpha,-lambda,-n_nonull_coef,
         -compartment) |>
  sapply(
    FUN = function(x) {
      nonzero <- sum(x != 0)
      nonzero_prop <- nonzero / length(x)
      x_nonzero <- x[x != 0]
      min <- min(x_nonzero)
      mean <- mean(x_nonzero)
      median <- median(x_nonzero)
      max <- max(x_nonzero)
      results <- c(
        "gene" = names(x),
        "nonzero_values" = nonzero,
        "nonzero_prop" = nonzero_prop,
        "min_coef" = min,
        "max_coef" = max,
        "mean_coef" = mean,
        "median_coef" = median
      )
      results
    }
  ) |>
  t()

# Gene names
gene_names <- repcv_coef_summary_tumor |>
  rownames()

# Replace with 0 values
# Add gene names
# Add compartmnent
# Add overall sign of association with the outcome
# Set order of columns
repcv_coef_summary_tumor <- repcv_coef_summary_tumor |>
  as_tibble() |>
  mutate_all(
    .funs = function(x) {
      replace(x, !is.finite(x), 0)
    }
  ) |>
  mutate(
    gene = gene_names,
    compartment = "tumor",
    overall_sign = case_when(
      max_coef <= 0 ~ "negative association with LTR",
      min_coef >= 0 ~ "positive association with LTR",
      .default = "non-conclusive"
    )
  ) |>
  select(
    gene,
    compartment,
    nonzero_values,
    nonzero_prop,
    min_coef,
    max_coef,
    mean_coef,
    median_coef,
    overall_sign
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot",
    "model_tumor",
    "repcv_pred_tumor",
    "repcv_coef_tumor",
    "repcv_coef_summary_tumor"
  )
))

# 2. STROMAL COMPARTMENT ----

# Load results
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/4_OSIRESP_long_term_response_biomarkers_stroma.RData"
)

# Save object with results using a different name
model_stroma <- model_results

# Prediction results
repcv_pred_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$prediccions
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "stroma")

# Coefficients per gene
repcv_coef_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$coefficients
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "stroma")

# Check dimensions. Total number of elastic net models is 100 rep x 5 fold = 500
repcv_coef_stroma |> dim()

# Results per gene:
# - number of models where selected
# - proportion of the total of 500 models
# - minimum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - maximum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - mean coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - median coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
repcv_coef_summary_stroma <- repcv_coef_stroma |>
  select(-`(Intercept)`,-rep,-testfold,-alpha,-lambda,-n_nonull_coef,
         -compartment) |>
  sapply(
    FUN = function(x) {
      nonzero <- sum(x != 0)
      nonzero_prop <- nonzero / length(x)
      x_nonzero <- x[x != 0]
      min <- min(x_nonzero)
      mean <- mean(x_nonzero)
      median <- median(x_nonzero)
      max <- max(x_nonzero)
      results <- c(
        "gene" = names(x),
        "nonzero_values" = nonzero,
        "nonzero_prop" = nonzero_prop,
        "min_coef" = min,
        "max_coef" = max,
        "mean_coef" = mean,
        "median_coef" = median
      )
      results
    }
  ) |>
  t()

# Gene names
gene_names <- repcv_coef_summary_stroma |>
  rownames()

# Replace with 0 values
# Add gene names
# Add compartmnent
# Add overall sign of association with the outcome
# Set order of columns
repcv_coef_summary_stroma <- repcv_coef_summary_stroma |>
  as_tibble() |>
  mutate_all(
    .funs = function(x) {
      replace(x, !is.finite(x), 0)
    }
  ) |>
  mutate(
    gene = gene_names,
    compartment = "stroma",
    overall_sign = case_when(
      max_coef <= 0 ~ "negative association with LTR",
      min_coef >= 0 ~ "positive association with LTR",
      .default = "non-conclusive"
    )
  ) |>
  select(
    gene,
    compartment,
    nonzero_values,
    nonzero_prop,
    min_coef,
    max_coef,
    mean_coef,
    median_coef,
    overall_sign
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot",
    "model_tumor",
    "repcv_pred_tumor",
    "repcv_coef_tumor",
    "repcv_coef_summary_tumor",
    "model_stroma",
    "repcv_pred_stroma",
    "repcv_coef_stroma",
    "repcv_coef_summary_stroma"
  )
))

# 3. JOINT ANALYSIS -----

# Load results
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/5_OSIRESP_long_term_response_biomarkers_joint.RData"
)

# Save object with results using a different name
model_tumor_stroma <- model_results

# Prediction results
repcv_pred_tumor_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$prediccions
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor_stroma")

# Coefficients per gene
repcv_coef_tumor_stroma <- lapply(
  1:length(repcv_results),
  FUN = function(rep) {
    lapply(
      1:nfolds,
      FUN = function(fold) {
        repcv_results[[rep]][[fold]]$coefficients
      }
    ) |>
      bind_rows()
  }
) |>
  bind_rows() |>
  mutate(compartment = "tumor_stroma")

# Check dimensions. Total number of elastic net models is 100 rep x 5 fold = 500
repcv_coef_tumor_stroma |> dim()

# Results per gene:
# - number of models where selected
# - proportion of the total of 500 models
# - minimum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - maximum coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - mean coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
# - median coefficient value among those that are not zero (i.e., among those
# from models in which the gene is selected)
repcv_coef_summary_tumor_stroma <- repcv_coef_tumor_stroma |>
  select(-`(Intercept)`,-rep,-testfold,-alpha,-lambda,-n_nonull_coef,
         -compartment) |>
  sapply(
    FUN = function(x) {
      nonzero <- sum(x != 0)
      nonzero_prop <- nonzero / length(x)
      x_nonzero <- x[x != 0]
      min <- min(x_nonzero)
      mean <- mean(x_nonzero)
      median <- median(x_nonzero)
      max <- max(x_nonzero)
      results <- c(
        "gene" = names(x),
        "nonzero_values" = nonzero,
        "nonzero_prop" = nonzero_prop,
        "min_coef" = min,
        "max_coef" = max,
        "mean_coef" = mean,
        "median_coef" = median
      )
      results
    }
  ) |>
  t()

# Gene names
gene_names <- repcv_coef_summary_tumor_stroma |>
  rownames()

# Replace with 0 values
# Add gene names
# Add compartmnent
# Add overall sign of association with the outcome
# Set order of columns
repcv_coef_summary_tumor_stroma <- repcv_coef_summary_tumor_stroma |>
  as_tibble() |>
  mutate_all(
    .funs = function(x) {
      replace(x, !is.finite(x), 0)
    }
  ) |>
  mutate(
    gene = gene_names,
    compartment = "tumor_stroma",
    overall_sign = case_when(
      max_coef <= 0 ~ "negative association with LTR",
      min_coef >= 0 ~ "positive association with LTR",
      .default = "non-conclusive"
    )
  ) |>
  select(
    gene,
    compartment,
    nonzero_values,
    nonzero_prop,
    min_coef,
    max_coef,
    mean_coef,
    median_coef,
    overall_sign
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot",
    "model_tumor",
    "repcv_pred_tumor",
    "repcv_coef_tumor",
    "repcv_coef_summary_tumor",
    "model_stroma",
    "repcv_pred_stroma",
    "repcv_coef_stroma",
    "repcv_coef_summary_stroma",
    "model_tumor_stroma",
    "repcv_pred_tumor_stroma",
    "repcv_coef_tumor_stroma",
    "repcv_coef_summary_tumor_stroma"
  )
))

# 4. BIND RESULTS ----

# For graphical representation
# - prediction results
repcv_pred <- list(repcv_pred_tumor, repcv_pred_stroma, repcv_pred_tumor_stroma) |>
  bind_rows() |>
  mutate(compartment = compartment |>
           factor(
             levels = c("tumor", "stroma", "tumor_stroma"),
             labels = c("tumor", "stroma", "tumor and stroma"),
             ordered = TRUE
           ))

# - coefficients
repcv_coef <- list(repcv_coef_tumor, repcv_coef_stroma) |>
  bind_rows() |>
  mutate(compartment = compartment |>
           factor(
             levels = c("tumor", "stroma", "tumor_stroma"),
             labels = c("tumor", "stroma", "tumor and stroma"),
             
             ordered = TRUE
           ))

# - coefficients summary
repcv_coef_summary <- list(repcv_coef_summary_tumor, repcv_coef_summary_stroma) |>
  bind_rows() |>
  mutate(
    gene = str_remove_all(gene, "`"),
    compartment = compartment |>
      factor(
        levels = c("tumor", "stroma", "tumor_stroma"),
        labels = c("tumor", "stroma", "tumor and stroma"),
        ordered = TRUE
      ),
    overall_sign = overall_sign |>
      factor(
        levels = c(
          "negative association with LTR",
          "positive association with LTR",
          "non-conclusive"
        ),
        ordered = TRUE
      )
  )

# 5. FIGURES ----

# Lollipop plot: robust biomarkers ----

gene_inclusion_plot_lr <- ggplot(repcv_coef_summary |>
                                   filter(compartment != "tumor and stroma" &
                                            nonzero_prop >= 0.40)) +
  geom_segment(
    aes(
      x = reorder_within(gene, nonzero_prop, compartment),
      xend = reorder_within(gene, nonzero_prop, compartment),
      y = nonzero_prop,
      yend = 0,
      color = overall_sign
    ),
    size = 0.05
  ) +
  geom_point(aes(
    x = reorder_within(gene, nonzero_prop, compartment),
    y = nonzero_prop,
    color = overall_sign
  ),
  size = 1) +
  scale_x_reordered() +
  scale_color_manual(
    values = c(
      "negative association with LTR" = "brown1",
      "positive association with LTR" = "darkviolet"
    ),
    labels = c("negative association with LTR", "positive association with LTR")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent
  ) +
  labs(x = "", y = "selection rate", color = "") +
  coord_flip() +
  facet_wrap(~ compartment, scales = "free") +
  # guides(color = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 6),
    axis.text = element_text(color = "black", size = 5),
    axis.ticks.y = element_line(size = 0.1),
    axis.ticks.x = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    legend.spacing.y = unit(0, "lines"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0),
    strip.background = element_rect(size = 0.25)
  )

# Ridge plot: coefficients/effects across all repetitions ----

# Genes to subset
# (those with a percentage of inclusion >= 40%)
genes_to_subset_tumor <- repcv_coef_summary_tumor |>
  filter(nonzero_prop >= 0.40) |>
  select(gene) |>
  pull()

genes_to_subset_tumor |>
  length()

genes_to_subset_stroma <- repcv_coef_summary_stroma |>
  filter(nonzero_prop >= 0.40) |>
  select(gene) |>
  pull()

genes_to_subset_stroma |>
  length()

# Data in long format
# Filter selected genes based on percentage of inclusion
# Filter null values for regression coefficients
# Remove "`" from gene names (for plot)

# Check colnames (for pivot_longer)
repcv_coef |> colnames() |> head(n = 20)
repcv_coef |> colnames() |> tail(n = 20)

repcv_coef_long <- repcv_coef |>
  pivot_longer(`(Intercept)`:TNFSF12,
               names_to = "gene",
               values_to = "coef") |>
  filter((compartment == "tumor" &
            gene %in% genes_to_subset_tumor) |
           (compartment == "stroma" &
              gene %in% genes_to_subset_stroma)
  ) |>
  filter(coef != 0) |>
  mutate(gene = str_remove_all(gene, "`"))

# Add mean coefficient per gene
repcv_coef_long <- repcv_coef_long |>
  left_join(repcv_coef_summary |>
              select(gene, compartment, mean_coef),
            by = join_by(gene, compartment))


# Ridge plot of non-zero regression coefficient values per gen
gene_coef_plot_lr <- ggplot(repcv_coef_long, aes(
  x = coef,
  y = reorder(gene, mean_coef),
  fill = mean_coef
)) +
  geom_vline(xintercept = 0, lwd = 0.05) +
  geom_density_ridges(
    linewidth = 0.05,
    color = "black",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0, height = 0),
    point_shape = "|",
    point_size = 0.5,
    alpha = 0.7
  ) +
  geom_density_ridges(linewidth = 0.05, color = "black") +
  scale_fill_gradientn(
    colours = c("brown1", "white", "darkviolet"),
    values = scales::rescale(c(-1, 0, 1)),
    limits = c(-1.5, 1.5),
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  ) +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4, by = 1)) +
  facet_wrap(~ compartment, scales = "free") +
  labs(x = "regression coefficient", y = "", fill = "mean value") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(hjust = 0.5),
    line = element_line(linewidth = 0.05, color = "black"),
    text = element_text(size = 6),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    strip.background = element_rect(color = "black", linewidth = 0.25),
    legend.key.height = unit(0.01, "npc")
  )

# AUC values ----

repcv_auc <- repcv_pred |>
  select(compartment, rep, test_fold, auc_testfold) |>
  unique()

# Check dimensions
repcv_auc |>
  group_by(compartment) |>
  summarise(n = n())

alpha <- 0.05

# Get mean values per repetition
repcv_auc_mean <- repcv_auc |>
  group_by(compartment, rep) |>
  summarise(auc_testfold = mean(auc_testfold))

# Check dimensions
repcv_auc_mean |>
  group_by(compartment) |>
  summarise(n = n())

# Get point estimate and confidence intervals
repcv_auc_stats <- repcv_auc_mean |>
  group_by(compartment) |>
  summarise(
    mean = mean(auc_testfold),
    ci_low = quantile(auc_testfold, probs = alpha / 2),
    ci_upp = quantile(auc_testfold, probs = 1 - alpha / 2),
    label = paste0(
      format(round(mean, digits = 2), nsmall = 2),
      " (",
      format(round(ci_low, digits = 2), nsmall = 2),
      ", ",
      format(round(ci_upp, digits = 2), nsmall = 2),
      ")"
    )
  )

repcv_auc_stats

dodge_width <- 0.1

auc_plot_lr <- ggplot(data = repcv_auc_mean, aes(x = compartment, y = auc_testfold, fill = compartment)) +
  stat_halfeye(point_interval = NULL, alpha = 0.8) +
  geom_jitter(
    data = repcv_auc_mean,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(compartment) - dodge_width,
      y = auc_testfold,
      color = compartment
    ),
    height = 0,
    width = 0.025,
    alpha = 0.5,
    size = 0.25
  ) +
  geom_point(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(x = compartment, y = mean),
    size = 1
  ) +
  geom_errorbar(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(
      x = compartment,
      y = mean,
      ymin = ci_low,
      ymax = ci_upp
    ),
    width = 0.1,
    linewidth = 0.25
  ) +
  scale_y_continuous(limits = c(0.5, 0.9),
                     breaks = seq(0.5, 0.9, by = 0.05)) +
  scale_x_discrete(expand = c(0.25, 0)) +
  # scale_color_manual(values = c("tumor" = "#08CC3C",
  #                               "stroma" = "#121DEE",
  #                               "tumor and stroma" = "#B2001A")) +
  # scale_fill_manual(values = c("tumor" = "#08CC3C",
  #                              "stroma" = "#121DEE",
  #                              "tumor and stroma" = "#B2001A")) +
  # scale_fill_manual(values = c(
  #   "tumor" = "#2EAD73",
  #   "stroma" = "#121DEE",
  #   "tumor and stroma" = "#B2001A"
  # )) +
  # scale_color_manual(values = c(
  #   "tumor" = "#2EAD73",
  #   "stroma" = "#121DEE",
  #   "tumor and stroma" = "#B2001A"
  # )) +
  scale_color_manual(values = c(
    "tumor" = "#66DEAA",
    "stroma" = "#4A5FFF",
    "tumor and stroma" = "#993399"
  )) +
  scale_fill_manual(values = c(
    "tumor" = "#66DEAA",
    "stroma" = "#4A5FFF",
    "tumor and stroma" = "#993399"
  )) +
  geom_label(
    data = repcv_auc_stats,
    inherit.aes = FALSE,
    aes(x = compartment, y = mean, label = label),
    hjust = -0.15,
    # vjust = 0.2,
    show.legend = FALSE,
    size = 1.5,
    label.size = 0.1
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0),
    text = element_text(size = 6),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25)
  ) +
  labs(x = "", y = "AUC estimates from repeated cross-validation") +
  guides(color = guide_legend(nrow = 2))

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot",
    "gene_inclusion_plot_lr",
    "auc_plot_lr",
    "gene_coef_plot_lr"
  )
))

# GENE SIGNATURE OF LTR ----

# Load gene signature (RF prediction model) and validation results

# Spatially informed RF prediction model
spatial_rf_list <- readRDS(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/6_OSIRESP_long_term_response_signature_spatial.rds")
spatial_rf <- spatial_rf_list$spatial_rf

# Validation results for both the spatially resolved and the pseudobulk gene
# signature
load(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/7_OSIREAL_long_term_response_signature_validation.RData")

# Partial Dependence Plots (PDPs) ----

pdp <- lapply(
  spatial_rf$xvar.names,
  FUN = function(var) {
    plot.variable(
      spatial_rf,
      xvar.names = var,
      partial = TRUE,
      oob = TRUE,
      npts = 50,
      target = "Yes"
    )
  }
)

pdp_custom <- lapply(
  1:length(pdp),
  FUN = function(ind) {
    # Extract plot data for the gene
    plot_data_df <- tibble(x = pdp[[ind]]$pData[[1]]$x) |>
      left_join(
        tibble(
          x = pdp[[ind]]$pData[[1]]$x.uniq,
          yhat = pdp[[ind]]$pData[[1]]$yhat,
          yhat_se = pdp[[ind]]$pData[[1]]$yhat.se
        )
      )
    
    
    ceiling_dec <- function(x, level = 1) {
      power <- 10 ^ level
      ceiling(x * power) / power
    }
    
    floor_dec <- function(x, level = 1) {
      power <- 10 ^ level
      floor(x * power) / power
    }
    
    xmin <- plot_data_df$x |> min() |> floor_dec()
    xmax <- plot_data_df$x |> max() |> ceiling_dec()
    
    # Custom ggplot with similar style to plot.variable, but controlling axis limits
    ggplot(plot_data_df, aes(x = x, y = yhat)) +
      geom_point(color = "red", size = 0.1) +
      geom_line(
        data = plot_data_df |>
          filter(!is.na(yhat)),
        inherit.aes = FALSE,
        aes(x = x, y = yhat),
        color = "black",
        linetype = "dashed",
        linewidth = 0.05
      ) +
      geom_line(
        data = plot_data_df |>
          filter(!is.na(yhat_se)),
        inherit.aes = FALSE,
        aes(x = x, y = yhat + yhat_se),
        color = "red",
        linetype = "dotted",
        linewidth = 0.05
      ) +
      geom_line(
        data = plot_data_df |>
          filter(!is.na(yhat_se)),
        inherit.aes = FALSE,
        aes(x = x, y = yhat - yhat_se),
        color = "red",
        linetype = "dotted",
        linewidth = 0.05
      ) +
      geom_rug(sides = "b", linewidth = 0.2) +
      theme_classic() +
      theme(text = element_text(size = 6),
            line = element_line(linewidth = 0.075, color = "black")) +
      labs(x = sub("_", " in ", pdp[[ind]]$pData[[1]]$xvar.names),
           y = "LTR probability") +
      scale_y_continuous(limits = c(0.25, 0.50)) +
      scale_x_continuous(
        limits = c(xmin, xmax),
        breaks = seq(xmin, xmax, length.out = 5),
        labels = round(seq(xmin, xmax, length.out = 5), digits = 1)
      )
  }
)

# Spearman's correlation from PDPs ----

spearman_corr <- lapply(
  1:length(pdp),
  FUN = function(ind) {
    x <- pdp[[ind]]$pData[[1]]$x.uniq
    y <- pdp[[ind]]$pData[[1]]$yhat
    
    estimate <- cor(x, y, method = "spearman")
    names(estimate) <- pdp[[ind]]$pData[[1]]$xvar.names
    
    tibble(
      gene = sub("_", " in ", names(estimate)),
      spearman_corr = estimate,
      sign = estimate |> sign() |>
        factor(
          levels = c(-1, 1),
          labels = c(
            "negative association with LTR",
            "positive association with LTR"
          )
        )
    )
  }
) |>
  bind_rows()

# Plot
spearman_pdp_plot <- ggplot(spearman_corr, aes(
  x = spearman_corr,
  y = reorder(gene, spearman_corr),
  color = sign
)) +
  geom_vline(xintercept = 0, lwd = 0.05) +
  geom_segment(aes(
    x = spearman_corr,
    xend = 0,
    y = reorder(gene, spearman_corr),
    yend = gene,
    color = sign
  ),
  size = 0.05,
  ) +
  geom_point(aes(x = spearman_corr, y = gene, color = sign, ), size = 0.5, ) +
  scale_color_manual(
    values = c(
      "positive association with LTR" = "darkviolet",
      "negative association with LTR" = "brown1"
    )
  ) +
  labs(
    x = expression("spearman's " * rho * " from partial dependence plot"),
    y = "",
    color = "",
    fill = ""
  ) +
  scale_size(limits = c(0, 1)) +
  theme_classic() +
  # guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 6),
    line = element_line(linewidth = 0.075, color = "black"),
    legend.key.height = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-1, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2))

# ROC curves ----

# 1. Spatially resolved gene signature ----

oob_roc <- roc(
  response = spatial_rf_list$spatial_rf_oob$long_response,
  predictor = spatial_rf_list$spatial_rf_oob$oob_pred
)

seed <- 123456
n.boot <- 10000

# Spatially resolved prediction model
set.seed(seed)
oob_auc_ci <- ci.auc(oob_roc, method = "boot", boot.n = n.boot)

roc_labs <- tibble(
  validation = c("OSIRESP cohort (out-of-bag)", "OSIREAL cohort"),
  auc = format(
    c(oob_roc$auc, advanced_auc$auc),
    digits = 2,
    nsmall = 2
  ),
  ci_low = format(
    c(oob_auc_ci[1], advanced_auc_ci[1]),
    digits = 2,
    nsmall = 2
  ),
  ci_upp = format(
    c(oob_auc_ci[3], advanced_auc_ci[3]),
    digits = 2,
    nsmall = 2
  ),
  label = paste0("AUC = ", auc, " (", ci_low, ", ", ci_upp, ")")
)

roc_plot <- ggroc(list(
  "OSIRESP cohort (out-of-bag)" = oob_roc,
  "OSIREAL cohort" = advanced_auc
),
size = 0.25) +
  geom_segment(aes(
    x = 1,
    y = 0,
    xend = 0,
    yend = 1
  ),
  color = "gray37",
  linewidth = 0.05) +
  coord_fixed() +
  geom_label(
    inherit.aes = FALSE,
    data = roc_labs,
    aes(
      x = c(0.95, 0.57),
      y = c(0.98, 0.825),
      color = validation,
      label = label
    ),
    hjust = 0,
    size = 2,
    label.size = 0.1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "OSIREAL cohort" = "#EF453E",
    "OSIRESP cohort (out-of-bag)" = "red3"
  )) +
  labs(x = "specificity", y = "sensitivity", color = "") +
  # guides(label = "none",
  #        color = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 6),
    line = element_line(linewidth = 0.075, color = "black"),
    legend.key.height = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-1, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2))

# 2. Spatially resolved + pseudobulk gene signatures ----

roc_supp_labs <- tibble(
  validation = c(
    "Spatially resolved gene signature (OSIREAL cohort)",
    "Pseudobulk gene signature (OSIREAL cohort)"
  ),
  auc = format(
    c(advanced_auc$auc, advanced_auc_pseudo$auc),
    digits = 2,
    nsmall = 2
  ),
  ci_low = format(
    c(advanced_auc_ci[1], advanced_auc_ci_pseudo[1]),
    digits = 2,
    nsmall = 2
  ),
  ci_upp = format(
    c(advanced_auc_ci[3], advanced_auc_ci_pseudo[3]),
    digits = 2,
    nsmall = 2
  ),
  label = paste0("AUC = ", auc, " (", ci_low, ", ", ci_upp, ")")
)

roc_pseudo_plot <- ggroc(
  list(
    "Spatially resolved gene signature (OSIREAL cohort)" = advanced_auc,
    "Pseudobulk gene signature (OSIREAL cohort)" = advanced_auc_pseudo
  ),
  size = 0.25
) +
  geom_segment(aes(
    x = 1,
    y = 0,
    xend = 0,
    yend = 1
  ),
  color = "gray37",
  linewidth = 0.05) +
  coord_fixed() +
  geom_label(
    inherit.aes = FALSE,
    data = roc_supp_labs,
    aes(
      x = c(0.75, 0.57),
      y = c(0.95, 0.825),
      color = validation,
      label = label
    ),
    hjust = 0,
    size = 2,
    label.size = 0.1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Pseudobulk gene signature (OSIREAL cohort)" = "#EF453E",
      "Spatially resolved gene signature (OSIREAL cohort)" = "red3"
    )
  ) +
  labs(x = "specificity", y = "sensitivity", color = "") +
  # guides(label = "none",
  #        color = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 6),
    line = element_line(linewidth = 0.075, color = "black"),
    legend.key.height = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-1, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2))

# Kaplan-Meier curves ----

# 1. OSIRESP cohort ----
# (OOB predictions)
osiresp_pred_oob <- spatial_rf_list$spatial_rf_oob
osiresp_pred_oob_q2 <- quantile(osiresp_pred_oob$oob_pred, probs = 0.5)

osiresp_pred_oob <- osiresp_pred_oob |>
  mutate(
    pred_q2 = ifelse(
      oob_pred >= osiresp_pred_oob_q2,
      "LTR prediction \u2265 median",
      "LTR prediction < median"
    ) |>
      factor(
        levels = c("LTR prediction \u2265 median", "LTR prediction < median")
      )
  )

osiresp_sf <- survfit(
  Surv(pfs_time_months, disease_progression_to_osimertinib) ~ pred_q2,
  data = osiresp_pred_oob |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.type = "log-log"
)

osiresp_sf_med <- quantile(osiresp_sf, probs = 0.5)$quantile[, "50"] |>
  round(digits = 2)

osiresp_km <- ggsurvplot(
  osiresp_sf,
  data = osiresp_pred_oob |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "",
  legend.labs = levels(osiresp_pred_oob$pred_q2),
  # legend.labs = paste0("C", 1:4),
  palette = c("darkviolet", "brown1"),
  break.x.by = 12,
  xlim = c(0, 72),
  ggtheme = theme_classic(base_size = 6),
  risk.table.fontsize = 2,
  pval.size = 2,
  size = 0.25,
  censor.size = 2,
  xlab = "Time (months)",
  ylab = "Progression-free survival",
  tables.y.text = FALSE
)

osiresp_km$plot <- osiresp_km$plot +
  geom_segment(
    data = data.frame(
      x = osiresp_sf_med,
      xend = osiresp_sf_med,
      y = rep(0, 2),
      yend = rep(0.5, 2)
    ),
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    size = 0.1,
    linetype = "dashed"
  ) +
  geom_segment(
    data = data.frame(
      x = 0,
      xend = max(osiresp_sf_med),
      y = 0.5,
      yend = 0.5
    ),
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    size = 0.1,
    linetype = "dashed"
  ) +
  theme(
    axis.line = element_line(size = 0.1),
    axis.ticks = element_line(size = 0.1),
    text = element_text(size = 6)
  )

osiresp_km$table <- osiresp_km$table +
  theme(
    plot.title = element_text(size = 6),
    axis.line = element_line(size = 0.1),
    axis.ticks = element_line(size = 0.1),
    text = element_text(size = 6)
  )


# 2. OSIREAL cohort ----
osireal_pred_q2 <- quantile(advanced_pred$pred, probs = 0.5)

advanced_pred <- advanced_pred |>
  mutate(
    pred_q2 = ifelse(
      pred >= osireal_pred_q2,
      "LTR prediction \u2265 median",
      "LTR prediction < median"
    ) |>
      factor(
        levels = c("LTR prediction \u2265 median", "LTR prediction < median")
      )
  )

# Global
osireal_sf <- survfit(
  Surv(pfs_time_months, disease_progression_to_osimertinib) ~ pred_q2,
  data = advanced_pred |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.type = "log-log"
)

osireal_sf_med <- quantile(osireal_sf, probs = 0.5)$quantile[, "50"] |>
  round(digits = 2)

osireal_km <- ggsurvplot(
  osireal_sf,
  data = advanced_pred |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "",
  legend.labs = levels(advanced_pred$pred_q2),
  # legend.labs = paste0("C", 1:4),
  palette = c("darkviolet", "brown1"),
  break.x.by = 12,
  xlim = c(0, 72),
  ggtheme = theme_classic(base_size = 6),
  risk.table.fontsize = 2,
  pval.size = 2,
  size = 0.25,
  censor.size = 2,
  xlab = "Time (months)",
  ylab = "Progression-free survival",
  tables.y.text = FALSE
)

osireal_km$plot <- osireal_km$plot +
  geom_segment(
    data = data.frame(
      x = osireal_sf_med,
      xend = osireal_sf_med,
      y = rep(0, 2),
      yend = rep(0.5, 2)
    ),
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    size = 0.1,
    linetype = "dashed"
  ) +
  geom_segment(
    data = data.frame(
      x = 0,
      xend = max(osireal_sf_med),
      y = 0.5,
      yend = 0.5
    ),
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    size = 0.1,
    linetype = "dashed"
  ) +
  theme(
    axis.line = element_line(size = 0.1),
    axis.ticks = element_line(size = 0.1),
    text = element_text(size = 6)
  )

osireal_km$table <- osireal_km$table +
  theme(
    plot.title = element_text(size = 6),
    axis.line = element_line(size = 0.1),
    axis.ticks = element_line(size = 0.1),
    text = element_text(size = 6)
  )

# Remove objects from work space
rm(list = setdiff(
  ls(),
  c(
    "dodge_width",
    "alpha",
    "gene_inclusion_plot",
    "gene_coef_plot",
    "auc_plot",
    "gene_inclusion_plot_lr",
    "gene_coef_plot_lr",
    "auc_plot_lr",
    "pdp_custom",
    "spearman_pdp_plot",
    "osiresp_km",
    "osireal_km",
    "roc_plot",
    "roc_pseudo_plot"
  )
))

# FIGURES ----

# Figure 2 ----

# Diagram
diagram <- readPNG("./Results/Figures/Biorender_diagrams/Final_figures/Methods_biomarkers.png")

diagram_ggplot <- ggdraw() +
  draw_image(diagram)

dsp_fig <- plot_grid(
  NULL,
  diagram_ggplot,
  NULL,
  nrow = 1,
  rel_widths = c(0.025, 1, 0.1),
  labels = c("A", "", "")
)

# Results from SIS + ENET analyses (Cox regression)
enet_fig <- plot_grid(
  gene_inclusion_plot,
  auc_plot,
  nrow = 1,
  labels = c("B", "C"),
  rel_widths = c(1.6, 1)
)

# Results from SIS + ENET analyses (logistic regression)
auc_fig <- plot_grid(NULL,
                     auc_plot_lr,
                     nrow = 1,
                     rel_widths = c(0.02, 1))

# Results from SIS + ENET analyses (logistic regression) +
# spearman correlation coefficient estimated from Partial Depedence Plots
# capturing direction of association with LTR for each gene comprised in the
# spatially resolved signature of LTR
lr_pdp <- plot_grid(auc_fig,
                    spearman_pdp_plot,
                    nrow = 1,
                    labels = c("D", "E"))

# Partial Dependence Plots from Random Forest model for genes comprised in the
# spatially resolved signature of LTR
# (selected genes for Figure 2)
pdp_subset <- plot_grid(
  pdp_custom[[4]] +
    scale_y_continuous(limits = c(0.25, 0.46)) +
    scale_x_continuous(limits = c(5, 8.5), breaks = seq(5, 8.5, by = 0.5)),
  pdp_custom[[1]] +
    scale_y_continuous(limits = c(0.25, 0.46)) +
    scale_x_continuous(limits = c(5, 8.5), breaks = seq(5, 8.5, by = 0.5)),
  pdp_custom[[9]] +
    scale_y_continuous(limits = c(0.25, 0.46)) +
    scale_x_continuous(limits = c(5, 8.5), breaks = seq(5, 8.5, by = 0.5)),
  pdp_custom[[7]] +
    scale_y_continuous(limits = c(0.25, 0.46)) +
    scale_x_continuous(limits = c(5, 8.5), breaks = seq(5, 8.5, by = 0.5)),
  nrow = 2
)

# Arrange previous plots
lr_pdp_spearman <- plot_grid(
  lr_pdp,
  pdp_subset,
  nrow = 1,
  labels = c("", "F"),
  rel_widths = c(1.6, 1)
)

# Kaplan-Meier plots
osiresp_km_grob <- arrangeGrob(osiresp_km$plot,
                               osiresp_km$table,
                               ncol = 1,
                               heights = c(2, 0.75))

osireal_km_grob <- arrangeGrob(osireal_km$plot,
                               osireal_km$table,
                               ncol = 1,
                               heights = c(2, 0.75))

# ROC curves
roc_fig <- plot_grid(NULL,
                     roc_plot,
                     nrow = 1,
                     rel_widths = c(0.015, 1))

# Arrange previous plots
roc_km_fig <- plot_grid(
  roc_fig,
  osiresp_km_grob,
  osireal_km_grob,
  nrow = 1,
  labels = c("G", "H", "I")
)

# Figure 2
fig_2 <- plot_grid(
  enet_fig,
  lr_pdp_spearman,
  roc_km_fig,
  nrow = 3,
  rel_heights = c(0.8, 0.8, 1)
)

# Figure 2 (adding diagram)
fig_2_diagram <- plot_grid(dsp_fig,
                           fig_2,
                           nrow = 2,
                           rel_heights = c(0.25, 1)) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save figure
ggsave(
  "./Results/Figures/Final_figures/Figure_2.png",
  plot = fig_2_diagram,
  units = "cm",
  width = 21,
  height = 29,
  dpi = 600
)

# Supp Figure S1 ----

# Arrange previous plots
fig_s1 <- plot_grid(
  gene_coef_plot,
  gene_inclusion_plot_lr,
  gene_coef_plot_lr,
  nrow = 3,
  labels = c("A", "B", "C")
)

# Save figure
ggsave(
  "./Results/Figures/Final_figures/Figure_S1.png",
  plot = fig_s1,
  units = "cm",
  width = 21,
  height = 29,
  dpi = 600
)

# Supp Figure S2 ----

# Partial Dependence Plots for all genes in the spatially resolved signature of
# LTR

# Margin adjustment
pdp_custom_axis <- lapply(
  pdp_custom,
  FUN = function(plot) {
    plot <- plot +
      scale_x_continuous(limits = c(5, 11.5), breaks = seq(5, 11, by = 1)) +
      scale_y_continuous(limits = c(0.25, 0.50))
    
    plot
    
  }
)

# Plots
fig_s2_pdp <- plot_grid(plotlist = pdp_custom_axis, ncol = 4) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Add left margin to figure
fig_s2_pdp_margin <- plot_grid(NULL,
                               fig_s2_pdp,
                               nrow = 1,
                               rel_widths = c(0.025, 1))

# ROC curve (spatially and pseudo-bulk signature)
fig_s2_roc_bulk <- plot_grid(roc_pseudo_plot, NULL, nrow = 1)

fig_s2 <- plot_grid(
  fig_s2_pdp_margin,
  fig_s2_roc_bulk,
  nrow = 2,
  labels = c("A", "B"),
  rel_heights = c(2.75 / 3, 1 / 2)
) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(
  "./Results/Figures/Final_figures/Figure_S2.png",
  plot = fig_s2,
  units = "cm",
  width = 21,
  height = 26,
  dpi = 600
)
