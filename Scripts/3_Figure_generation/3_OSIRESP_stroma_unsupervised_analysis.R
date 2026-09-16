# DESCRIPTION ----
# Figures 3 and S3

# LIBRARIES ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggdist)
library(pheatmap)
library(survival)
library(survminer)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggplotify)
library(rstatix)
library(png)
library(ggtext)

# COLORS ----

# Color for expression levels in heatmap
heat_col <- colorRampPalette(c("blue3", "white", "red3"))(100)

clust_col <- c(
  "C1_A" = "#6666E6",
  "C2_A" = "darkorchid",
  "C3_A" = "chocolate2",
  "C4_A" = "brown3"
)

# Annotation colors for heatmap
ann_heat_col <- list(
  cluster = clust_col,
  stage = c(
    "II" = "#D0E4FF",
    "III" = "royalblue",
    "IV" = "#0000CD"
  ),
  `tobbaco history` = c(
    "Never smoker" = "#FFD6D6",
    "Former smoker" = "indianred",
    "Active smoker" = "#C82333"
  ),
  `EGFR mutation` = c(
    "Exon 19 deletion" = "#D3A84E",
    "Exon 21 L858R" = "lightcoral",
    "Other" = "royalblue"
  ),
  `TP53 mutation` = c(
    "Yes" = "darkred",
    "No" = "darkgoldenrod",
    "Non-available" = "lightgray"
  )
)

# LOAD EXPRESSION DATA ----

expr_data <- readRDS("./Data/OSIRESP_cohort/Processed/OSIRESP_vst_medianvalues.rds")

# LOAD RESULTS ----

load(file = "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/8_OSIRESP_stroma_unsupervised_analysis.RData")

# UMAP (stroma) ----

# Cluster size and percentage
cluster_size <- leid_res$cluster |>
  table()

cluster_names <- names(cluster_size)

cluster_perc <- round(cluster_size / sum(cluster_size) * 100, digits = 2)

cluster_size <- cluster_size |>
  as.vector()

cluster_perc <- cluster_perc |>
  as.vector()

# Modify cluster names in leiden clustering results
leid_res <- leid_res |>
  mutate(cluster = paste0(cluster, "_A"))

# Set shapes for UMAP plot
leid_shapes <- c(22:25)
names(leid_shapes) <- names(clust_col)

# UMAP colored by cluster
umap_leid <- ggplot(leid_res,
                    aes(
                      x = umap_1,
                      y = umap_2,
                      color = cluster,
                      fill = cluster,
                      shape = cluster
                    )) +
  geom_point(size = 0.5) +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  scale_shape_manual(values = leid_shapes) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-4, 4)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 8),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing = unit(0, "pt"),
    legend.box.margin = margin(-10, 0, -5, 0)
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")

# HEATMAP GSVA ----

# 1. Stroma compartment ----

# Enrichment scores variance per gene set
gsva_stroma_var <- apply(gsva_stroma, MARGIN = 1, FUN = var)

# First quartile for variance values
gsva_stroma_var_p <- quantile(gsva_stroma_var, prob = 0.35)

# Seleccion of top 75% variable gene sets
gsva_stroma_var <- gsva_stroma_var[gsva_stroma_var > gsva_stroma_var_p]

# Sample annotation
annot_stroma <- data.frame(
  patient_id = leid_res$patient_id,
  cluster = leid_res$cluster,
  `tobbaco history` = leid_res$tobacco_history,
  `EGFR mutation` = leid_res$egfr_mutation_type,
  `TP53 mutation` = leid_res$tp53_mut,
  stage = leid_res$stage,
  check.names = FALSE
)

# Sample annotation: order by cluster
annot_stroma <- annot_stroma |>
  arrange(cluster)

# Columns order
annot_stroma <- annot_stroma |>
  select(patient_id,
         cluster,
         `tobbaco history`,
         `EGFR mutation`,
         `TP53 mutation`,
         stage)

# Sample annotation: set rownames
rownames(annot_stroma) <- annot_stroma$patient_id

# Color legend
breaks_gsva_stroma <- seq(-0.7, 0.7, length.out = 100)

# Gap cols
gap_cols_stroma <- head(as.numeric(cumsum(table(
  annot_stroma$cluster
))), -1)

# Heatmap of GSVA scores (all gene sets)
gsva_stroma_heatmap <- pheatmap(
  gsva_stroma[, rownames(annot_stroma)],
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_method = "ward.D",
  annotation_col = annot_stroma[, 2:5],
  show_rownames = TRUE,
  show_colnames = FALSE,
  breaks = breaks_gsva_stroma,
  color = heat_col,
  annotation_colors = ann_heat_col,
  gaps_col = gap_cols_stroma,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 4,
  fontsize_row = 3,
  cellwidth = 2.25,
  cellheight = 2.5,
  border_color = NA
)$gtable |>
  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Heatmap including gene sets with more score variability
gsva_stroma_heatmap_subset <- pheatmap(
  gsva_stroma[names(gsva_stroma_var), rownames(annot_stroma)],
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_method = "ward.D",
  annotation_col = annot_stroma[, 2:5],
  show_rownames = TRUE,
  show_colnames = FALSE,
  breaks = breaks_gsva_stroma,
  color = heat_col,
  annotation_colors = ann_heat_col,
  gaps_col = gap_cols_stroma,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 4,
  fontsize_row = 4,
  cellwidth = 1.35,
  cellheight = 3.5,
  border_color = NA
)$gtable |>
  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# 2. Tumor compartment ----

annot_tumor <- annot_stroma |>
  filter(patient_id %in% colnames(gsva_tumor))

annot_stroma |>
  dim()

annot_tumor |>
  dim()

annot_tumor$cluster |> table()

# Gap cols
gap_cols_tumor <- head(as.numeric(cumsum(table(
  annot_tumor$cluster
))), -1)

# Heatmap of GSVA scores (all gene sets)
gsva_heat_tumor <- pheatmap(
  gsva_tumor[, rownames(annot_tumor)],
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_method = "ward.D",
  annotation_col = annot_stroma[, -1],
  show_rownames = TRUE,
  show_colnames = FALSE,
  breaks = breaks_gsva_stroma,
  color = heat_col,
  annotation_colors = ann_heat_col,
  gaps_col = gap_cols_tumor,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 4,
  fontsize_row = 3,
  cellwidth = 2.25,
  cellheight = 2.5,
  border_color = NA
)$gtable |>
  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# GSVA SCORES FOR GENE SETS OF INTEREST ----

# Tumor compartment ----

gsva_tumor_clust <- gsva_tumor |>
  t() |>
  as.tibble() |>
  mutate(patient_id = colnames(gsva_tumor)) |>
  pivot_longer(cols = -patient_id,
               names_to = "gene_set",
               values_to = "gsva_score") |>
  mutate(compartment = "Tumor") |>
  inner_join(annot_tumor |>
               select(patient_id, cluster))

gs_tumor <- c(
  "PDGF Signaling",
  "Cordenonsi YAP Conserved Signature",
  "EGFR Signaling",
  "MET Signaling"
)

gsva_tumor_scores_plot <- ggplot(
  gsva_tumor_clust |>
    filter(gene_set %in% gs_tumor) |>
    mutate(gene_set = gene_set |>
             factor(levels = gs_tumor, ordered = TRUE)),
  aes(y = cluster, x = gsva_score, fill = cluster)
) +
  geom_violin(trim = FALSE,
              color = "white",
              alpha = 0.65) +
  scale_fill_manual(values = clust_col) +
  scale_color_manual(values = clust_col) +
  geom_boxplot(alpha = 0, size = 0.15) +
  geom_jitter(
    shape = 16,
    size = 0.25,
    color = "black",
    height = 0.05,
    width = 0
  ) +
  geom_vline(xintercept = 0,
             linetype = "dotted",
             lwd = 0.2) +
  labs(y = "", x = "GSVA enrichment score in tumor") +
  scale_x_continuous(limits = c(-1, 1)) +
  facet_wrap(~ gene_set) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.x = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0)
  )

# TOP MARKERS C4 A (stroma) ----

# Filter significantly DE genes with log2FC > 0 (that is, upregulated in
# cluster C4_A vs all others)
mrkr_tumor_c4_up <- mrkr_res_stroma$`C4_A cluster` |>
  filter(adj.P.Val < 0.05 & logFC > 0)

# Filter top 20 DE genes with log2FC > 0, based on log2FC magnitude
mrkr_tumor_c4_up_top <- mrkr_tumor_c4_up |>
  arrange(desc(logFC)) |>
  slice(1:50) |>
  mutate(contrast = "c4_vs_other") |>
  rownames_to_column(var = "gene")

# Color palette for plots
pal <- viridisLite::magma(256)
subset_pal <- pal[60:150]

# Top markers - C4_A subtype
markr_c4_up_plot <- ggplot(mrkr_tumor_c4_up_top |>
                             slice(1:20),
                           aes(
                             x = rev(reorder(gene, logFC)),
                             y = logFC,
                             # fill = adj.P.Val,
                             color = adj.P.Val
                           )) +
  geom_hline(
    yintercept = 1,
    color = "black",
    linetype = "dotted",
    lwd = 0.2
  ) +
  geom_segment(aes(
    x = rev(reorder(gene, logFC)),
    xend = rev(reorder(gene, logFC)),
    y = logFC,
    yend = 0,
    color = adj.P.Val,
  ),
  size = 0.15) +
  geom_point(aes(
    x = rev(reorder(gene, logFC)),
    y = logFC,
    color = adj.P.Val,
    size = logFC
  ), ) +
  geom_point(show.legend = FALSE) +
  scale_color_continuous(palette = subset_pal) +
  geom_richtext(
    aes(label = gene),
    angle = 45,
    fill = "white",
    size = 1.25,
    label.size = 0.1,
    vjust = -0.5,
    hjust = 0,
    label.padding = unit(0.1, "lines"),
    color = "black"
  ) +
  facet_wrap(
    ~ contrast,
    scales = "free_y",
    labeller = labeller(contrast = c(c4_vs_other = "cluster C4_A")),
    ncol = 4
  ) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5)) +
  scale_color_gradientn(
    colours = subset_pal,
    limits = c(0, 0.05),
    breaks = seq(0.01, 0.05, 0.01)
  ) +
  scale_size_continuous(range = c(1, 3)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.1),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0)
  ) +
  guides(color = guide_colorbar(barwidth = 9, barheight = 0.3), size = "none") +
  labs(x = "",
       y = expression("log"[2] * "FC in stroma"),
       color = "adj. p")

# MCAF GENES (Stroma) ----

# Join with cluster results
leid_expr <- leid_res |>
  left_join(expr_data) |>
  mutate(compartment = compartment |>
           factor(
             levels = c("Tumor", "Stroma"),
             labels = c("Tumor", "Stroma"),
             ordered = TRUE
           ))

mcaf_genes <- c(
  "COL11A1",
  "PDPN",
  "MMP14",
  "COL1A1",
  "VIM",
  "SMA",
  "ACTA2",
  "S100A4",
  "DCN",
  "LUM",
  "VCAN",
  "COL14A1",
  "FBLN1",
  "FBLN2",
  "SMOC",
  "LOX",
  "LOXL1",
  "CXCL14",
  "COL6A3",
  "FN1",
  "COL3A1",
  "SPON2",
  "COL5A1",
  "PDGFRA",
  "INHBA"
)

mcaf_genes |> length()

# mCAF markers comprised in the CTA assay
mcaf_genes <- intersect(mcaf_genes, colnames(expr_data))

# Data for plots
mcaf_genes_stroma <- leid_expr |>
  filter(compartment == "Stroma") |>
  select(c("cluster", mcaf_genes))

# Long format
mcaf_genes_stroma_long <- mcaf_genes_stroma |>
  pivot_longer(cols = -cluster,
               names_to = "gene",
               values_to = "expression")

# p-value position (Y axis)
mcaf_y_pos <- mcaf_genes_stroma_long |>
  group_by(gene) |>
  summarise(y.position = max(expression, na.rm = TRUE) - 0.5)

# p-value formatting
mcaf_pvals <- mrkr_res_stroma$`C4_A cluster`[mcaf_genes, ] |>
  rownames_to_column(var = "gene") |>
  mutate(
    group1 = "C1_A",
    group2 = "C4_A",
    signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  inner_join(mcaf_y_pos)

# Plots
dodge_width <- 0.1

mcaf_boxplots <- ggplot(
  mcaf_genes_stroma_long |>
    mutate(
      cluster = recode(
        cluster,
        "C1" = "C1_A",
        "C2" = "C2_A",
        "C3" = "C3_A",
        "C4" = "C4_A"
      )
    ),
  aes(
    x = cluster,
    y = expression,
    color = cluster,
    fill = cluster
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.65) +
  geom_boxplot(
    width = 0.5,
    size = 0.1,
    alpha = 0,
    color = "black"
  ) +
  facet_wrap(~ gene, scales = "free") +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0)
  ) +
  stat_pvalue_manual(
    mcaf_pvals,
    label = "signif",
    tip.length = 0,
    bracket.size = 0.1,
    size = 2
  ) +
  labs(x = "", y = expression("log"[2] * " expression in stroma"))

# LOXL2 ----
dodge_width <- 0.1

# 1. Stroma compartment ----

# Data for plot
loxl2_stroma <- leid_expr |>
  filter(compartment == "Stroma") |>
  select(c("cluster", LOXL2))

# Long format
loxl2_stroma_long <- loxl2_stroma |>
  pivot_longer(cols = -cluster,
               names_to = "gene",
               values_to = "expression")

# p-value position (Y axis)
loxl2_y_pos <- loxl2_stroma_long |>
  group_by(gene) |>
  summarise(y.position = max(expression, na.rm = TRUE) - 0.1)

# p-value formatting
loxl2_pvals <- mrkr_res_stroma$`C4_A cluster`["LOXL2", ] |>
  rownames_to_column(var = "gene") |>
  mutate(
    group1 = "C1_A",
    group2 = "C4_A",
    signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  inner_join(loxl2_y_pos)

# Plot
dodge_width <- 0.1

loxl2_stroma <- ggplot(
  loxl2_stroma_long |>
    mutate(
      cluster = recode(
        cluster,
        "C1" = "C1_A",
        "C2" = "C2_A",
        "C3" = "C3_A",
        "C4" = "C4_A"
      )
    ),
  aes(
    x = cluster,
    y = expression,
    color = cluster,
    fill = cluster
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.65) +
  geom_jitter(
    data = loxl2_stroma_long,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(as.factor(cluster)) - dodge_width,
      y = expression,
      color = cluster
    ),
    height = 0,
    width = 0.01,
    size = 0.01
  ) +
  geom_boxplot(
    width = 0.5,
    size = 0.1,
    alpha = 0,
    color = "black"
  ) +
  facet_wrap(~ gene, scales = "free") +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0)
  ) +
  stat_pvalue_manual(
    loxl2_pvals,
    label = "signif",
    tip.length = 0,
    bracket.size = 0.1,
    size = 2
  ) +
  labs(x = "", y = expression("log"[2] * " expression in stroma"))

# 2. Tumor compartment ----

# Data for plot
loxl2_tumor <- leid_expr |>
  filter(compartment == "Tumor") |>
  select(c("cluster", LOXL2))

# Long format
loxl2_tumor_long <- loxl2_tumor |>
  pivot_longer(cols = -cluster,
               names_to = "gene",
               values_to = "expression")

# p-value position
loxl2_y_pos <- loxl2_tumor_long |>
  group_by(gene) |>
  summarise(y.position = max(expression, na.rm = TRUE) - 0.1)

# p-value formatting
loxl2_pvals <- mrkr_res_stroma$`C4_A cluster`["LOXL2", ] |>
  rownames_to_column(var = "gene") |>
  mutate(
    group1 = "C1_A",
    group2 = "C4_A",
    signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  inner_join(loxl2_y_pos)

# Plot
dodge_width <- 0.1

loxl2_tumor <- ggplot(
  loxl2_tumor_long |>
    mutate(
      cluster = recode(
        cluster,
        "C1" = "C1_A",
        "C2" = "C2_A",
        "C3" = "C3_A",
        "C4" = "C4_A"
      )
    ),
  aes(
    x = cluster,
    y = expression,
    color = cluster,
    fill = cluster
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.65) +
  geom_jitter(
    data = loxl2_tumor_long,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(as.factor(cluster)) - dodge_width,
      y = expression,
      color = cluster
    ),
    height = 0,
    width = 0.01,
    size = 0.01
  ) +
  geom_boxplot(
    width = 0.5,
    size = 0.1,
    alpha = 0,
    color = "black"
  ) +
  facet_wrap(~ gene, scales = "free") +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0)
  ) +
  stat_pvalue_manual(
    loxl2_pvals,
    label = "signif",
    tip.length = 0,
    bracket.size = 0.1,
    size = 2
  ) +
  labs(x = "", y = expression("log"[2] * " expression in tumor"))

# T AND B CELLS ----

# 1. Stroma compartment ----

# Selected markers
t_and_b_cells_genes <- c("CD8A", "CD8B", "CD19", "CD38")

# Data for plots
t_and_b_cells_genes_stroma <- leid_expr |>
  filter(compartment == "Stroma") |>
  select(c("cluster", t_and_b_cells_genes))

# Long format
t_and_b_cells_genes_stroma_long <- t_and_b_cells_genes_stroma |>
  pivot_longer(cols = -cluster,
               names_to = "gene",
               values_to = "expression") |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# p-value position
t_and_b_cells_y_pos <- t_and_b_cells_genes_stroma_long |>
  group_by(gene) |>
  summarise(y.position = max(expression, na.rm = TRUE) - 0.1) |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# p-value formatting
t_and_b_cells_pvals <- mrkr_res_stroma$`C4_A cluster`[t_and_b_cells_genes, ] |>
  rownames_to_column(var = "gene") |>
  mutate(
    group1 = "C1_A",
    group2 = "C4_A",
    signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  inner_join(t_and_b_cells_y_pos) |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# Plots
dodge_width <- 0.1

t_and_b_cells_stroma <- ggplot(
  t_and_b_cells_genes_stroma_long |>
    mutate(
      cluster = recode(
        cluster,
        "C1" = "C1_A",
        "C2" = "C2_A",
        "C3" = "C3_A",
        "C4" = "C4_A"
      )
    ),
  aes(
    x = cluster,
    y = expression,
    color = cluster,
    fill = cluster
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.65) +
  geom_jitter(
    data = t_and_b_cells_genes_stroma_long,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(as.factor(cluster)) - dodge_width,
      y = expression,
      color = cluster
    ),
    height = 0,
    width = 0.01,
    size = 0.01
  ) +
  geom_boxplot(
    width = 0.5,
    size = 0.1,
    alpha = 0,
    color = "black"
  ) +
  facet_wrap(~ gene, scales = "free", nrow = 1) +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0)
  ) +
  stat_pvalue_manual(
    t_and_b_cells_pvals,
    label = "signif",
    tip.length = 0,
    bracket.size = 0.1,
    size = 2
  ) +
  labs(x = "", y = expression("log"[2] * " expression in stroma"))

# 2. Tumor compartment ----

# Data for plots
t_and_b_cells_genes_tumor <- leid_expr |>
  filter(compartment == "Tumor") |>
  select(c("cluster", t_and_b_cells_genes))

# Long format
t_and_b_cells_genes_tumor_long <- t_and_b_cells_genes_tumor |>
  pivot_longer(cols = -cluster,
               names_to = "gene",
               values_to = "expression") |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# p-value position (Y axis)
t_and_b_cells_y_pos <- t_and_b_cells_genes_tumor_long |>
  group_by(gene) |>
  summarise(y.position = max(expression, na.rm = TRUE) - 0.1) |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# p-value formatting
t_and_b_cells_pvals <- mrkr_res_tumor$`C4_A cluster`[t_and_b_cells_genes, ] |>
  rownames_to_column(var = "gene") |>
  mutate(
    group1 = "C1_A",
    group2 = "C4_A",
    signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  inner_join(t_and_b_cells_y_pos) |>
  mutate(gene = gene |> factor(levels = t_and_b_cells_genes, ordered = TRUE))

# Plots
dodge_width <- 0.1

t_and_b_cells_tumor <- ggplot(
  t_and_b_cells_genes_tumor_long |>
    mutate(
      cluster = recode(
        cluster,
        "C1" = "C1_A",
        "C2" = "C2_A",
        "C3" = "C3_A",
        "C4" = "C4_A"
      )
    ),
  aes(
    x = cluster,
    y = expression,
    color = cluster,
    fill = cluster
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.65) +
  geom_jitter(
    data = t_and_b_cells_genes_tumor_long,
    inherit.aes = FALSE,
    aes(
      x = as.numeric(as.factor(cluster)) - dodge_width,
      y = expression,
      color = cluster
    ),
    height = 0,
    width = 0.01,
    size = 0.01
  ) +
  geom_boxplot(
    width = 0.5,
    size = 0.1,
    alpha = 0,
    color = "black"
  ) +
  facet_wrap(~ gene, scales = "free", nrow = 1) +
  scale_color_manual(values = clust_col) +
  scale_fill_manual(values = clust_col) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0)
  ) +
  stat_pvalue_manual(
    t_and_b_cells_pvals,
    label = "signif",
    tip.length = 0,
    bracket.size = 0.1,
    size = 2
  ) +
  labs(x = "", y = expression("log"[2] * " expression in tumor"))

# KAPLAN-MEIER ----

leid_sf <- survfit(
  Surv(pfs_time_months, disease_progression_to_osimertinib) ~ cluster,
  data = leid_res |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.type = "log-log"
)

leid_sf

# Median PFS times
leid_sf_med <- quantile(leid_sf, probs = 0.5)$quantile[, "50"] |>
  round(digits = 2)

# Kaplan-Meier plot
leid_km <- ggsurvplot(
  leid_sf,
  data = leid_res |>
    mutate(
      disease_progression_to_osimertinib = ifelse(disease_progression_to_osimertinib == "Yes", 1, 0)
    ),
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "",
  legend.labs = paste0("C", 1:4, "_A"),
  palette = unname(clust_col),
  break.x.by = 12,
  xlim = c(0, 72),
  ggtheme = theme_classic(base_size = 8),
  risk.table.fontsize = 2,
  pval.size = 2,
  size = 0.25,
  censor.size = 2,
  xlab = "Time (months)",
  ylab = "Progression-Free Survival",
  tables.y.text = FALSE
)

# Change graphical parameters used in plot
leid_km$plot <- leid_km$plot +
  geom_segment(
    data = data.frame(
      x = leid_sf_med,
      xend = leid_sf_med,
      y = rep(0, 4),
      yend = rep(0.5, 4)
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
      xend = max(leid_sf_med),
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
    text = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, -10, 0)
  )

# Change graphical parameters used in plot (risk table)
leid_km$table <- leid_km$table +
  labs(x = "") +
  theme(
    plot.title = element_text(size = 8),
    axis.line = element_line(size = 0.1),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    text = element_text(size = 6)
  )

# Rearrange plot and risk table
leid_km <- arrangeGrob(leid_km$plot,
                       leid_km$table,
                       ncol = 1,
                       heights = c(3, 0.8))

# FIGURES ----

# Figure 3 ----

# Biorender diagram
diagram <- readPNG("./Results/Figures/Biorender_diagrams/Clustering_advanced.png")

diagram_ggplot <- ggdraw() +
  draw_image(diagram)

# UMAP and Kaplan-Meier curves
fig_umap_km <- plot_grid(umap_leid,
                         leid_km,
                         nrow = 1,
                         labels = c("B", "C"))

# Heatmap of GSVA scores (stroma compartment, most variable gene sets)
fig_gsva_heat <- plot_grid(NULL,
                           gsva_stroma_heatmap_subset,
                           nrow = 1,
                           rel_widths = c(0.075, 1))

# C4_A markers
fig_mrkr_c4 <- plot_grid(NULL,
                         markr_c4_up_plot,
                         nrow = 1,
                         rel_widths = c(0.025, 1))

# Arrange GSVA heatmap and C4_A markers
fig_gsva_stroma_c4_mrkr <- plot_grid(fig_gsva_heat,
                                     fig_mrkr_c4,
                                     nrow = 1,
                                     labels = c("D", "E"))
# Heatmap of GSVA scores (tumor comparment)
fig_gsva_scores <- plot_grid(NULL,
                             gsva_tumor_scores_plot,
                             nrow = 1,
                             rel_widths = c(0.025, 1))

# Arrange mCAF plots (stromal compartment) and GSVA scores (tumor compartment)
fig_mcaf_gsva_tumor <- plot_grid(mcaf_boxplots,
                                 fig_gsva_scores,
                                 nrow = 1,
                                 labels = c("F", "G"))

# Figure
fig <- plot_grid(fig_umap_km,
                 fig_gsva_stroma_c4_mrkr,
                 fig_mcaf_gsva_tumor,
                 nrow = 3)

# Figure, adding diagram
fig_diagram <- plot_grid(
  diagram_ggplot,
  fig,
  nrow = 2,
  labels = c("A", ""),
  rel_heights = c(0.4, 3)
) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(
  "./Results/Figures/Final_figures/Figure_3.png",
  plot = fig_diagram,
  units = "cm",
  width = 21,
  height = 29,
  dpi = 600
)

# Figure S3 ----

# Heatmap of GSVA scores (stroma, all gene sets)
fig_heat_stroma <- plot_grid(NULL,
                             gsva_stroma_heatmap,
                             nrow = 1,
                             rel_widths = c(0.06, 1))

# Arrange heatmap of GSVA scores and LOXL2 plot (stromal compartment)
fig_stroma <- plot_grid(
  fig_heat_stroma,
  loxl2_stroma,
  nrow = 1,
  rel_widths = c(1.75 / 3, 1.05 / 3),
  labels = c("A", "B")
)

# Heatmap of GSVA scores (tumor compartment, all gene sets)
fig_heat_tumor <- plot_grid(NULL,
                            gsva_heat_tumor,
                            nrow = 1,
                            rel_widths = c(0.06, 1))

# Arrange heatmap of GSVA scores and LOXL2 plot (tumor compartment)
fig_tumor <- plot_grid(
  fig_heat_tumor,
  loxl2_tumor,
  nrow = 1,
  rel_widths = c(1.75 / 3, 1.05 / 3),
  labels = c("D", "E")
)

# Figure
fig_s3 <- plot_grid(
  fig_stroma,
  t_and_b_cells_stroma,
  fig_tumor,
  t_and_b_cells_tumor,
  nrow = 4,
  rel_heights = c(1, 0.65, 1, 0.65),
  labels = c("", "C", "", "F")
) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(
  "./Results/Figures/Final_figures/Figure_S3.png",
  plot = fig_s3,
  units = "cm",
  width = 21,
  height = 29,
  dpi = 600
)
