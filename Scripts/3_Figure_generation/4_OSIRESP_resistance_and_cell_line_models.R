# DESCRIPTION:
# Figures 4, S4, and S5

# LIBRARIES ----
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggtext)
library(tidyr)
library(ggdist)
library(tidytext)
library(patchwork)
library(png)
library(ggplotify)
library(cowplot)
library(enrichplot)

# PATIENT SAMPLES ----

# Post-progression vs. pre-treatment paired patient tumor samples
load(
  "./Results/Intermediate/OSIRESP_and_OSIREAL_cohorts/9_OSIRESP_RB_postprogression_vs_pretreatment.RData"
)

# Color for expression levels in heatmap
col <- colorRampPalette(c("blue2", "white", "brown2"))(100)

# 1. Genes with higher logFC ----

fit_stroma_res <- fit_stroma_res |>
  rownames_to_column(var = "gene")

top_logFC_stroma <- fit_stroma_res |>
  filter(logFC > 0) |>
  arrange(desc(t)) |>
  slice(1:25)

# Color palette for plots
pal <- viridisLite::magma(256)
subset_pal <- pal[45:150]

# Top markers - all stromal subtypes
top_logFC_stroma_plot <- ggplot(top_logFC_stroma |>
                                  slice(1:20),
                                aes(
                                  x = reorder(gene, logFC, decreasing = TRUE),
                                  y = logFC,
                                  color = P.Value
                                )) +
  geom_hline(
    yintercept = 1,
    color = "black",
    linetype = "dotted",
    lwd = 0.2
  ) +
  geom_segment(aes(
    x = reorder(gene, logFC, decreasing = TRUE),
    xend = reorder(gene, logFC, decreasing = TRUE),
    y = logFC,
    yend = 0,
    color = P.Value,
  ),
  size = 0.15) +
  geom_point(aes(
    x = reorder(gene, logFC, decreasing = TRUE),
    y = logFC,
    color = P.Value,
    size = logFC
  ), ) +
  geom_richtext(
    inherit.aes = FALSE,
    data = top_logFC_stroma |> slice(1:20),
    aes(
      y = logFC,
      x = reorder(gene, logFC, decreasing = TRUE),
      label = gene
    ),
    fill = "white",
    color = "black",
    angle = 45,
    size = 1.25,
    label.size = 0.1,
    vjust = -0.5,
    hjust = 0,
    label.padding = unit(0.1, "lines")
  ) +
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
  labs(y = "",
       x = expression("log"[2] * "FC in stroma"),
       fill = "p-value")

# Spaghetti plot markers stroma compartment ----

# Selected genes for plot
gene_subset <- c(
  "COL5A1",
  "COL5A2",
  "LOXL2",
  "COMP",
  "VCAN",
  "MMP11",
  "NID2",
  "PDPN",
  "PLAU",
  "PLAUR",
  "SPP1",
  "PDGFA"
)

# Prepare data
expr_paired_stroma_subset <- expr_paired_stroma |>
  select(c("patient_id", "biopsy_type", gene_subset)) |>
  pivot_longer(COL5A1:PDGFA, names_to = "gene", values_to = "value") |>
  mutate(biopsy_type = biopsy_type |>
           factor(
             levels = c("pre", "post"),
             labels = c("pre-treatment", "at disease progression")
           ))

# Add lines between paired pre- and post- samples
expr_paired_stroma_subset <- expr_paired_stroma_subset |>
  group_by(gene, patient_id) |>
  mutate(line_col = if_else(value[biopsy_type == "at disease progression"] >
                              value[biopsy_type == "pre-treatment"], "increase", "decrease")) |>
  ungroup()

p_stroma <- ggplot(
  expr_paired_stroma_subset,
  aes(
    x = biopsy_type,
    y = value,
    color = biopsy_type,
    fill = biopsy_type
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.25) +
  geom_line(mapping = aes(group = patient_id, color = line_col),
            lwd = 0.2) +
  facet_wrap(~ gene) +
  scale_fill_manual(values = c('gray50', 'brown2')) +
  scale_color_manual(
    values = c(
      `pre-treatment` = 'gray50',
      `at disease progression` = 'brown2',
      `increase` = 'brown2',
      `decrease` = 'gray50'
    )
  ) +
  geom_point(size = 0.6) +
  facet_wrap(~ gene, scales = "free") +
  labs(
    x = "",
    y = expression("log"[2] * " expression in stroma"),
    color = "biopsy",
    fill = "biopsy"
  ) +
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
    legend.box.margin = margin(-10, 0, -10, 0)
  ) +
  guides(color = guide_colorbar(barwidth = 9, barheight = 0.3), size = "none")

# Spaghetti plot mCAF markers (stroma) ----

# mCAF markers
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

# mCAF markers comprised in the CTA assay
mcaf_genes <- intersect(mcaf_genes, colnames(expr_paired_stroma))

# Prepare data for plot
expr_paired_stroma_mcaf <- expr_paired_stroma |>
  select(c("patient_id", "biopsy_type", mcaf_genes)) |>
  pivot_longer(COL11A1:INHBA, names_to = "gene", values_to = "value") |>
  mutate(biopsy_type = biopsy_type |>
           factor(
             levels = c("pre", "post"),
             labels = c("pre-treatment", "at disease progression")
           ))

# Add lines between paired pre- and post-samples
expr_paired_stroma_mcaf <- expr_paired_stroma_mcaf |>
  group_by(gene, patient_id) |>
  mutate(line_col = if_else(value[biopsy_type == "at disease progression"] >
                              value[biopsy_type == "pre-treatment"], "increase", "decrease")) |>
  ungroup()

# Plot
p_mcaf <- ggplot(
  expr_paired_stroma_mcaf,
  aes(
    x = biopsy_type,
    y = value,
    color = biopsy_type,
    fill = biopsy_type
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.25) +
  geom_line(mapping = aes(group = patient_id, color = line_col),
            lwd = 0.2) +
  facet_wrap(~ gene) +
  scale_fill_manual(values = c('gray50', 'brown2')) +
  scale_color_manual(
    values = c(
      `pre-treatment` = 'gray50',
      `at disease progression` = 'brown2',
      `increase` = 'brown2',
      `decrease` = 'gray50'
    )
  ) +
  geom_point(size = 0.65) +
  facet_wrap(~ gene, scales = "free") +
  labs(
    x = "",
    y = expression("log"[2] * " expression in stroma"),
    color = "biopsy",
    fill = "biopsy"
  ) +
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
    legend.box.margin = margin(-10, 0, -10, 0)
  ) +
  guides(color = guide_colorbar(barwidth = 9, barheight = 0.3), size = "none")

# Spaghetti plot immune markers (stroma) ----

# Selected genes
immune_subset <- c(
  "CD3D",
  "CD3E",
  "CD247",
  "TRBC1/2",
  "ZAP70",
  "CD2",
  "CD79A",
  "CD27",
  "CD69",
  "CD38",
  "PDCD1",
  "POU2AF1"
)

# Prepare data for plot
expr_paired_stroma_immune_subset <- expr_paired_stroma |>
  select(c("patient_id", "biopsy_type", immune_subset)) |>
  pivot_longer(CD3D:POU2AF1, names_to = "gene", values_to = "value") |>
  mutate(biopsy_type = biopsy_type |>
           factor(
             levels = c("pre", "post"),
             labels = c("pre-treatment", "at disease progression")
           ))

# Add lines between paired pre- and post- samples
expr_paired_stroma_immune_subset <- expr_paired_stroma_immune_subset |>
  group_by(gene, patient_id) |>
  mutate(line_col = if_else(value[biopsy_type == "at disease progression"] >
                              value[biopsy_type == "pre-treatment"], "increase", "decrease")) |>
  ungroup()

# Plot
p_immune_subset <- ggplot(
  expr_paired_stroma_immune_subset,
  aes(
    x = biopsy_type,
    y = value,
    color = biopsy_type,
    fill = biopsy_type
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.25) +
  geom_line(mapping = aes(group = patient_id, color = line_col),
            lwd = 0.2) +
  facet_wrap(~ gene) +
  scale_fill_manual(values = c('gray50', 'brown2')) +
  scale_color_manual(
    values = c(
      `pre-treatment` = 'gray50',
      `at disease progression` = 'brown2',
      `increase` = 'brown2',
      `decrease` = 'gray50'
    )
  ) +
  geom_point(size = 0.6) +
  facet_wrap(~ gene, scales = "free", ncol = 4) +
  labs(
    x = "",
    y = expression("log"[2] * " expression in stroma"),
    color = "biopsy",
    fill = "biopsy"
  ) +
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
    legend.box.margin = margin(-10, 0, -10, 0)
  ) +
  guides(color = guide_colorbar(barwidth = 9, barheight = 0.3), size = "none")

# Spaghetti plot markers tumor compartment ----

# Selected genes
gene_subset <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2",
  "PLOD2",
  "VCAN",
  "THBS1",
  "ACTA2",
  "ITGA2",
  "LAMC2",
  "LAMB3"
)

# Prepare data for plot
expr_paired_tumor_subset <- expr_paired_tumor |>
  select(c("patient_id", "biopsy_type", gene_subset)) |>
  pivot_longer(COL1A1:LAMB3, names_to = "gene", values_to = "value") |>
  mutate(biopsy_type = biopsy_type |>
           factor(
             levels = c("pre", "post"),
             labels = c("pre-treatment", "at disease progression")
           ))

# Add lines between paired pre- and post samples
expr_paired_tumor_subset <- expr_paired_tumor_subset |>
  group_by(gene, patient_id) |>
  mutate(line_col = if_else(value[biopsy_type == "at disease progression"] >
                              value[biopsy_type == "pre-treatment"], "increase", "decrease")) |>
  ungroup()

p_tumor <- ggplot(
  expr_paired_tumor_subset,
  aes(
    x = biopsy_type,
    y = value,
    color = biopsy_type,
    fill = biopsy_type
  )
) +
  stat_halfeye(point_interval = NULL, alpha = 0.25) +
  geom_line(mapping = aes(group = patient_id, color = line_col),
            lwd = 0.2) +
  facet_wrap(~ gene) +
  scale_fill_manual(values = c('gray50', 'brown2')) +
  scale_color_manual(
    values = c(
      `pre-treatment` = 'gray50',
      `at disease progression` = 'brown2',
      `increase` = 'brown2',
      `decrease` = 'gray50'
    )
  ) +
  geom_point(size = 0.65) +
  facet_wrap(~ gene, scales = "free") +
  labs(
    x = "",
    y = expression("log"[2] * " expression in tumor"),
    color = "biopsy",
    fill = "biopsy"
  ) +
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
    legend.box.margin = margin(-10, 0, -10, 0)
  ) +
  guides(color = guide_colorbar(barwidth = 9, barheight = 0.3), size = "none")

# 2. GSEA ----

# All enriched gene sets ----

# Prepare data for plot
gsea_patients <- gsea_stroma@result |>
  mutate(compartment = "stroma") |>
  bind_rows(gsea_tumor@result |>
              mutate(compartment = "tumor"))

# Plot
gsea_patients_plot <- ggplot(
  gsea_patients |>
    filter(p.adjust < 0.05),
  aes(
    x = NES,
    y = reorder_within(ID, NES, compartment),
    size = -log10(p.adjust),
    color = NES
  )
) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  scale_y_reordered() +
  geom_vline(xintercept = 0,
             linewidth = 0.2,
             linetype = "dotted") +
  geom_point(show.legend = FALSE, alpha = 0.85) +
  scale_color_continuous(palette = col) +
  scale_size_continuous(range = c(1, 2.5)) +
  facet_wrap(~ compartment, scales = "free") +
  labs(x = "Normalized Enrichment Score", y = "") +
  theme_classic() +
  theme(
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-15, 0, 0, 0),
    text = element_text(size = 6),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    strip.background = element_rect(size = 0.25)
  )

# Selected enriched gene sets  ----

# Arrange data by significance in GSEA results
fit_stroma_res <- fit_stroma_res |>
  arrange(desc(t))

fit_tumor_res <- fit_tumor_res |>
  arrange(desc(t))

# Function for running enrichment score plot
gsea_plot_fct <- function(object,
                          geneSetID,
                          fit_res,
                          title = str_replace_all(geneSetID, "_", " "),
                          subtitle = NULL) {
  res <- object@result |>
    filter(ID == geneSetID)
  
  label <- sprintf("    p = %.2e\nadj. p = %.2e", res$pvalue, res$p.adjust)
  
  NES <- res$NES
  
  gs <- enrichplot:::gsInfo(object = object, geneSetID = geneSetID)
  
  offset <- 0.1 * diff(range(gs$runningScore))
  
  if (NES >= 0) {
    idx <- which.max(gs$runningScore)
    y_nes <- gs$runningScore[idx] - offset
    hjust_nes = -0.25
  } else {
    idx <- which.min(gs$runningScore)
    y_nes <- gs$runningScore[idx] + offset
    hjust_nes = 0.5
  }
  
  x_nes <- gs$x[idx]
  
  p_es <-
    ggplot(gs, aes(x, runningScore)) +
    geom_line(size = 0.8,
              colour = "green",
              lwd = 0.1) +
    geom_hline(yintercept = 0,
               linetype = "dotted",
               lwd = 0.2) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label,
      hjust = 1.15,
      vjust = 2.5,
      size = 1.25,
    ) +
    geom_label(
      data = tibble(
        x = x_nes,
        y = y_nes,
        label = sprintf("NES = %.2f", NES)
      ),
      inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      size = 1.25,
      label.size = 0.1,
      hjust = hjust_nes
    ) +
    labs(
      y = "Running Enrichment Score",
      x = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(size = 0.1),
      axis.ticks.y = element_line(size = 0.1),
      text = element_text(size = 5.5),
      title = element_text(size = 5.5)
    )
  
  hits <- subset(gs, position == 1)
  
  p_hits <-
    ggplot(hits, aes(x = x)) +
    geom_segment(aes(xend = x, y = 0, yend = 1), linewidth = 0.05) +
    theme_void()
  
  if (!is.null(fit_res)) {
    p_fc <-
      ggplot(fit_res, aes(
        x = 1:nrow(fit_res),
        y = 1,
        fill = logFC
      )) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradient2(
        low = "blue3",
        mid = "white",
        high = "red3",
        midpoint = 0,
        name = expression(log[2] * FC),
        limits = c(-1, 1),
        breaks = seq(-1, 1, by = 0.5),
        oob = scales::squish,
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        text = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        legend.key.height = unit(0.1, "cm")
      )
    
    gsea_plot <-
      p_es /
      p_hits /
      p_fc +
      plot_layout(heights = c(5, 0.8, 0.4))
    
  }
  
  else{
    gsea_plot <-
      p_es /
      p_hits +
      plot_layout(heights = c(5, 0.8, 0.4))
  }
  
  gsea_plot
}

# Selected gene sets (stromal compartment)
stroma_p1 <- gsea_plot_fct(gsea_stroma,
                           "Matrix Remodeling and Metastasis",
                           fit_stroma_res,
                           subtitle = "Stroma")

stroma_p2 <- gsea_plot_fct(gsea_stroma, "TGF-beta Signaling", fit_stroma_res, subtitle = "Stroma")

stroma_p3 <- gsea_plot_fct(gsea_stroma, "T cells", fit_stroma_res, subtitle = "Stroma")

stroma_p4 <- gsea_plot_fct(gsea_stroma, "B cells", fit_stroma_res, subtitle = "Stroma")

# Selected gene sets (tumor comparment)
tumor_p1 <- gsea_plot_fct(gsea_tumor, "PDGF Signaling", fit_tumor_res, subtitle = "Tumor")

tumor_p2 <- gsea_plot_fct(gsea_tumor, "TGF-beta Signaling", fit_tumor_res, subtitle = "Tumor")

tumor_p3 <- gsea_plot_fct(gsea_tumor,
                          "Cordenonsi YAP Conserved Signature",
                          fit_tumor_res,
                          subtitle = "Tumor")

tumor_p4 <- gsea_plot_fct(gsea_tumor, "MET Signaling", fit_tumor_res, subtitle = "Tumor")

# CELL LINE RESISTANCE MODEL ----

load(
  "./Results/Intermediate/Cell_line_models/Osimertinib_resistance/1_Osimertinib_resistance_cell_lines_GSEA.RData"
)

fdr <- 0.05

# Bind results for HCC827
gsea_hcc827_df <- lapply(
  1:length(gsea_hcc827),
  FUN = function(ind) {
    gsea_hcc827[[ind]]@result |>
      mutate(parental = names(gsea_hcc827)[ind])
  }
) |>
  bind_rows()

# Bind results for HCC827 and H1975
gsea_react_df <- gsea_hcc827_df |>
  mutate(
    cell_line = "HCC827",
    parental = recode(
      parental,
      MC1 = "monoclonal MC1",
      MC2 = "monoclonal MC2",
      MC3 = "monoclonal MC3"
    )
  ) |>
  rbind(gsea_h1975@result |>
          mutate(cell_line = "H1975", parental = "polyclonal")) |>
  mutate(parental = parental |>
           factor(
             levels = c(
               "polyclonal",
               "monoclonal MC1",
               "monoclonal MC2",
               "monoclonal MC3"
             ),
             ordered = TRUE
           ))

# Significant gene sets for all monoclonal parental cell lines in HCC827
# and for H1975. Only positively enriched gene sets.
gsea_react_signif_all <- gsea_react_df |>
  filter(p.adjust < fdr & NES > 0) |>
  group_by(ID) |>
  summarise(count = n()) |>
  filter(count == 4)

# Bind results for HCC827 and H1975
gsea_react_df <- gsea_react_df |>
  filter(ID %in% gsea_react_signif_all$ID)

# Format gene set names, add enrichment sign (for plot)
gsea_react_df <- gsea_react_df |>
  mutate(
    ID = str_remove_all(ID, "REACTOME") |>
      str_replace_all("_", " "),
    enrich_sign = ifelse(sign(NES) == 1, "up", "down") |>
      factor()
  )

# Plot
gsea_react_plot <- ggplot(
  gsea_react_df,
  aes(
    x = NES,
    y = reorder(ID, NES),
    shape = parental,
    color = cell_line,
    group = cell_line
  )
) +
  geom_errorbar(
    aes(xmin = 0, xmax = NES, color = cell_line),
    width = 0,
    position = position_dodge(width = 0.25),
    size = 0.25
  ) +
  geom_point(
    aes(
      x = NES,
      y = ID,
      color = cell_line,
      group = cell_line
    ),
    position = position_dodge(width = 0.25),
    size = 1
  ) +
  scale_color_manual(values = c("H1975"   = "#EF453E", "HCC827"  = "#D73027")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    shape = "",
    color = ""
  ) +
  guides(color = guide_legend(order = 1), shape = "none") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.margin     = margin(0, 0, 0, 0),
    legend.box = "vertical",
    legend.key.height = unit(0.5, "cm"),
    legend.spacing.y = unit(0.005, "cm"),
    text = element_text(size = 6),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    legend.box.margin = margin(-5, 0, -5, 0),
    axis.text.y = element_text(size = 3.5)
  )

# YAP/TAZ ----

# HCC827: monoclonal 1 parental cell line
log2fc_hcc827_mc1 <- log2fc_hcc827$MC1 |>
  rename(logFC = log2fc) |>
  arrange(desc(logFC))

res_hcc827_mc1_yap <- gsea_plot_fct(
  gsea_hcc827[[1]],
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  log2fc_hcc827_mc1,
  title = "Cordenonsi YAP Conserved Signature",
  subtitle = "HCC827 (MC1)"
)

# HCC827: monoclonal 2 parental cell line
log2fc_hcc827_mc2 <- log2fc_hcc827$MC2 |>
  rename(logFC = log2fc) |>
  arrange(desc(logFC))

res_hcc827_mc2_yap <- gsea_plot_fct(
  gsea_hcc827[[2]],
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  log2fc_hcc827_mc2,
  title = "Cordenonsi YAP Conserved Signature",
  subtitle = "HCC827 (MC2)"
)

# HCC827: monoclonal 3 parental cell line
log2fc_hcc827_mc3 <- log2fc_hcc827$MC3 |>
  rename(logFC = log2fc) |>
  arrange(desc(logFC))

res_hcc827_mc3_yap <- gsea_plot_fct(
  gsea_hcc827[[3]],
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  log2fc_hcc827_mc3,
  title = "Cordenonsi YAP Conserved Signature",
  subtitle = "HCC827 (MC3)"
)

# H1975: polyclonal parental cell line
log2fc_h1975 <- log2fc_h1975 |>
  rename(logFC = log2fc) |>
  arrange(desc(logFC))

res_h1975_yap <- gsea_plot_fct(
  gsea_h1975,
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  title = "Cordenonsi YAP Conserved Signature",
  log2fc_h1975,
  subtitle = "H1975"
)

# CELL LINE PERSISTENCE MODEL ----

load(
  "./Results/Intermediate/Cell_line_models/Persistence/1_GSE193258_cell_lines_GSEA.RData"
)

# Bind results for HCC827 and H1975
gsea_react_gse_df <- gsea_h1975@result |>
  mutate(cell_line = "H1975") |>
  rbind(gsea_hcc827@result |>
          mutate(cell_line = "HCC827")) |>
  rbind(gsea_hcc2935@result |>
          mutate(cell_line = "HCC2935")) |>
  mutate(cell_line = cell_line |>
           factor(
             levels = c("H1975", "HCC827", "HCC2935"),
             ordered = TRUE
           ))

# Significant gene sets for all monoclonal parental cell lines in HCC827
# and for H1975. Only positively enriched gene sets.
gsea_react_gse_signif_all <- gsea_react_gse_df |>
  filter(p.adjust < fdr & NES > 0) |>
  group_by(ID) |>
  summarise(count = n()) |>
  filter(count == 3)

# Bind results for HCC827 and H1975
gsea_react_gse_df <- gsea_react_gse_df |>
  filter(ID %in% gsea_react_gse_signif_all$ID)

# Format gene set names, add enrichment sign (for plot)
gsea_react_gse_df <- gsea_react_gse_df |>
  mutate(
    ID = str_remove_all(ID, "REACTOME") |>
      str_replace_all("_", " "),
    enrich_sign = ifelse(sign(NES) == 1, "up", "down") |>
      factor()
  )

# Recode REGULATION OF INSULIN LIKE GROWTH FACTOR IGF TRANSPORT AND UPTAKE BY INSULIN LIKE GROWTH FACTOR BINDING PROTEINS IGFBPS
gsea_react_gse_df <- gsea_react_gse_df |>
  mutate(
    ID = recode(
      ID,
      ` REGULATION OF INSULIN LIKE GROWTH FACTOR IGF TRANSPORT AND UPTAKE BY INSULIN LIKE GROWTH FACTOR BINDING PROTEINS IGFBPS` = " REGULATION OF INSULIN LIKE GROWTH FACTOR"
    )
  )

# Plot
gsea_react_gse_plot <- ggplot(gsea_react_gse_df,
                              aes(
                                x = NES,
                                y = reorder(ID, NES),
                                color = cell_line,
                                group = cell_line
                              )) +
  geom_errorbar(
    aes(xmin = 0, xmax = NES, color = cell_line),
    width = 0,
    position = position_dodge(width = 0.5),
    size = 0.25
  ) +
  geom_point(
    aes(
      x = NES,
      y = ID,
      color = cell_line,
      group = cell_line
    ),
    position = position_dodge(width = 0.5),
    size = 0.7,
    shape = 16
  ) +
  # scale_color_manual(values = c("H1975" = "deeppink3",
  #                               "HCC827" = "deeppink4",
  #                               "HCC2935" = "darkmagenta")) +
  scale_color_manual(values = c(
    "H1975"   = "#EF453E",
    # bright coral red
    "HCC827"  = "#D73027",
    # vivid red
    "HCC2935" = "#A50F15"   # deep crimson
  )) +
  scale_x_continuous(limits = c(0, 3)) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    shape = "",
    color = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.margin     = margin(0, 0, 0, 0),
    # legend.box = "vertical",
    legend.spacing.y = unit(0.005, "cm"),
    legend.key.height = unit(0.15, "cm"),
    text = element_text(size = 6),
    axis.ticks.x = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
    axis.line = element_line(size = 0.1),
    legend.box.margin = margin(-5, 0, -5, 0),
    axis.text.y = element_text(size = 3.5)
  )

# YAP/TAZ ----

# H1975
fit_results_h1975 <- fit_results_h1975 |>
  arrange(desc(t))

per_h1975_yap <- gsea_plot_fct(
  gsea_h1975,
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  fit_results_h1975,
  title = "Cordenonsi YAP Conserved Signature",
  
  subtitle = "H1975"
)

# HCC827
fit_results_hcc827 <- fit_results_hcc827 |>
  arrange(desc(t))

per_hcc827_yap <- gsea_plot_fct(
  gsea_hcc827,
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  fit_results_hcc827,
  title = "Cordenonsi YAP Conserved Signature",
  
  subtitle = "HCC827"
)

# HCC2935
fit_results_hcc2935 <- fit_results_hcc2935 |>
  arrange(desc(t))

per_hcc2935_yap <- gsea_plot_fct(
  gsea_hcc2935,
  "CORDENONSI_YAP_CONSERVED_SIGNATURE",
  fit_results_hcc2935,
  title = "Cordenonsi YAP Conserved Signature",
  
  subtitle = "HCC2935"
)

# FIGURES ----

# Figure 4 ----

# Patient biopsies
diagram_pat <- readPNG("./Results/Figures/Biorender_diagrams/Final_figures/Paired_biopsies.png")
diagram_pat <- ggdraw() +
  draw_image(diagram_pat)

# Resistance model
diagram_cell <- readPNG("./Results/Figures/Biorender_diagrams/Final_figures/Cell_lines_resistance.png")

diagram_cell <- ggdraw() +
  draw_image(diagram_cell)

# Running enrichment score plots (stromal compartment)
fig_stroma <- plot_grid(stroma_p1,
                        stroma_p2,
                        stroma_p3,
                        stroma_p4,
                        nrow = 2,
                        ncol = 2)

# Arrange mCAF marker plots and running enrichment score plots (stromal comparment)
fig_gene_stroma <- plot_grid(p_mcaf,
                             fig_stroma,
                             nrow = 1,
                             labels = c("B", "C"))

# Running enrichment score plots (tumor compartment)
fig_tumor <- plot_grid(tumor_p1,
                       tumor_p2,
                       tumor_p3,
                       tumor_p4,
                       nrow = 2,
                       ncol = 2)

# Arrange marker plots and running enrichment score plots (tumor compartment)
fig_gene_tumor <- plot_grid(p_tumor,
                            fig_tumor,
                            nrow = 1,
                            labels = c("D", "E"))

# Add left margin to diagram (osimertinib resistance cell line models)
fig_cell_diagram <- plot_grid(diagram_cell,
                              NULL,
                              nrow = 2,
                              rel_heights = c(1, 0.4))

# Arrange diagram and plot displaying enriched gene sets
fig_cell <- plot_grid(fig_cell_diagram,
                      gsea_react_plot,
                      nrow = 1,
                      labels = c("F", "G"))

# Add right margin to diagram (paired pre-treatment and post-progression tumor biopsies)
pat_fig <- plot_grid(diagram_pat,
                     NULL,
                     nrow = 1,
                     rel_widths = c(1, 0.2))

# Figure
fig_4 <- plot_grid(
  diagram_pat,
  fig_gene_stroma,
  NULL,
  fig_gene_tumor,
  NULL,
  fig_cell,
  nrow = 6,
  labels = c("A", "", "", "", "", ""),
  rel_heights = c(0.7, 1.5, 0.05, 1.5, 0.05, 1.25)
) |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(
  "./Results/Figures/Final_figures/Figure_4.png",
  plot = fig_4,
  units = "cm",
  width = 21,
  height = 26,
  dpi = 600
)

# Suplementary Figure S4 ----

# Selected top markers in stromal compartment
p_markers <- plot_grid(p_stroma,
                       p_immune_subset,
                       nrow = 1,
                       labels = c("A", "B"))

# Arrange marker plot and GSEA plot
fig_pat <- plot_grid(p_markers,
                     gsea_patients_plot,
                     nrow = 2,
                     labels = c("", "C"))

# YAP signature plots
fig_yap_res <- plot_grid(
  res_h1975_yap,
  res_hcc827_mc1_yap,
  res_hcc827_mc2_yap,
  res_hcc827_mc3_yap,
  nrow = 1,
  labels = c("D", "", "", "")
)

# Figure
fig_s4 <- plot_grid(fig_pat,
                    fig_yap_res,
                    # fig_yap_per,
                    nrow = 2,
                    rel_heights = c(3, 1))  |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(
  "./Results/Figures/Final_figures/Figure_S4.png",
  plot = fig_s4,
  units = "cm",
  width = 21,
  height = 21,
  dpi = 600
)

# Supplementary Figure S5 ----

# Diagram
diagram_cell_per <- readPNG(
  "./Results/Figures/Biorender_diagrams/Final_figures/Cell_lines_persistence.png"
)

diagram_cell_per <- ggdraw() +
  draw_image(diagram_cell_per)

# Add margins to diagram
fig_cell_diagram_per <- plot_grid(
  NULL,
  diagram_cell_per,
  NULL,
  NULL,
  nrow = 2,
  rel_heights = c(1, 0.55),
  rel_widths = c(0.025, 1)
)

# Arrange diagram and GSEA results
fig_cell_per <- plot_grid(
  fig_cell_diagram_per,
  gsea_react_gse_plot,
  nrow = 1,
  labels = c("A", "B")
)

# YAP signature plots
fig_yap_per <- plot_grid(
  per_h1975_yap,
  per_hcc827_yap,
  per_hcc2935_yap,
  nrow = 1,
  labels = c("C", "", "")
)

# Figure
fig_s5 <- plot_grid(fig_cell_per,
                    fig_yap_per,
                    nrow = 2,
                    rel_heights = c(1.5, 1))  |>
  as.ggplot() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(
  "./Results/Figures/Final_figures/Figure_S5.png",
  plot = fig_s5,
  units = "cm",
  width = 21,
  height = 12,
  dpi = 600
)
