# DESCRIPTION ----
# Figure 1


# LIBRARIES ----
library(png)
library(cowplot)
library(pdftools)
library(magick)
library(ggplot2)
library(ggplotify)

# BIORENDER DIAGRAMS ----

# Load diagrams
diagram_1 <- readPNG("./Results/Figures/Biorender_diagrams/Final_figures/Methods_1.png")
diagram_2 <- readPNG("./Results/Figures/Biorender_diagrams/Final_figures/Methods_2.png")

# Convert to cowplot drawable object
diagram_1_plot <- ggdraw() +
  draw_image(diagram_1)

diagram_2_plot <- ggdraw() +
  draw_image(diagram_2)

# GEOMX DSP IMAGES ----

# Load PNG files 
dsp_img_1 <- image_read("H12O_32_multiple_ROIs_1.png")
dsp_img_2 <- image_read("H12O_31_multiple_ROIs_1.png")
dsp_img_3 <- image_read("H12O_32_ROI004_nosegments_1.png")
dsp_img_4 <- image_read("H12O_32_ROI004_segments_1.png")

# Crop to same width/length
info1 <- image_info(dsp_img_1)
info2 <- image_info(dsp_img_2)
info3 <- image_info(dsp_img_3)
info4 <- image_info(dsp_img_4)

ratio1 <- info1$width / info1$height
ratio2 <- info2$width / info2$height
ratio3 <- info3$width / info3$height
ratio4 <- info4$width / info4$height

target_ratio <- min(ratio1, ratio2, ratio3, ratio4)

crop_to_ratio <- function(img, target_ratio) {
  info <- image_info(img)
  w <- info$width
  h <- info$height
  current_ratio <- w / h
  
  if (current_ratio > target_ratio) {
    # too wide: crop width (remove right side)
    new_w <- as.integer(h * target_ratio)
    # keep left edge fixed
    x_offset <- 0
    geometry <- paste0(new_w, "x", h, "+", x_offset, "+0")
    
  } else {
    # too tall: crop height (still center vertically)
    new_h <- as.integer(w / target_ratio)
    y_offset <- as.integer((h - new_h) / 2)
    geometry <- paste0(w, "x", new_h, "+0+", y_offset)
  }
  
  image_crop(img, geometry)
}

dsp_img_1_crop <- crop_to_ratio(dsp_img_1, target_ratio)
dsp_img_2_crop <- crop_to_ratio(dsp_img_2, target_ratio)
dsp_img_3_crop <- crop_to_ratio(dsp_img_3, target_ratio)
dsp_img_4_crop <- crop_to_ratio(dsp_img_4, target_ratio)

# Convert to cowplot drawable object
dsp_img_1_plot <- ggdraw() +
  draw_image(dsp_img_1_crop) +
  theme(plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))

dsp_img_2_plot <- ggdraw() +
  draw_image(dsp_img_2_crop) +
  theme(plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))

dsp_img_3_plot <- ggdraw() +
  draw_image(dsp_img_3_crop) +
  theme(plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))

dsp_img_4_plot <- ggdraw() +
  draw_image(dsp_img_4_crop) +
  theme(plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))

# Create legend for colors (panCK, SYTO12, and CD45)
biomarkers <- data.frame(x = 1:3,
                         y = 1:3,
                         marker = factor(
                           c("panCK", "SYTO13", "CD45"),
                           levels = c("panCK", "SYTO13", "CD45")
                         ))

biomarkers_plot <- ggplot(biomarkers, aes(
  x = x,
  y = y,
  color = marker,
  shape = marker
)) +
  geom_point(size = 0) +
  scale_color_manual(values = c(
    "panCK" = "#08CC3C",
    "SYTO13" = "#121DEE",
    "CD45" = "#D21E71"
  )) +
  # shape = 15 for solid square
  scale_shape_manual(values = rep(15, 3)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))

# Extract legend
biomarkers_legend <- get_legend(biomarkers_plot)
biomarkers_legend <- biomarkers_legend |>
  as.ggplot()

# Crop margins
biomarkers_legend <- cowplot::ggdraw(biomarkers_legend) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Add legend to GeoMx DSP images
dsp_img_1_plot_legend <- plot_grid(
  biomarkers_legend,
  dsp_img_1_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

dsp_img_2_plot_legend <- plot_grid(
  biomarkers_legend,
  dsp_img_2_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

dsp_img_3_plot_legend <- plot_grid(
  biomarkers_legend,
  dsp_img_3_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

dsp_img_4_plot_legend <- plot_grid(
  biomarkers_legend,
  dsp_img_4_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# FIGURE ----

# Arrange legends for GeoMx DSP images
fig_dsp_img_1 <- plot_grid(
  dsp_img_1_plot_legend,
  dsp_img_2_plot_legend,
  nrow = 1,
  labels = c("B", "C")
)

# Arrange GeoMx DSP images
fig_dsp_img_2 <- plot_grid(
  NULL,
  dsp_img_3_plot,
  NULL,
  dsp_img_4_plot,
  NULL,
  nrow = 1,
  rel_widths = c(0.5, 1, 0.025, 1, 0.5)
)

# Arrange GeoMx DSP images and legends
fig_dsp_img <- plot_grid(biomarkers_legend,
                         fig_dsp_img_2,
                         ncol = 1,
                         rel_heights = c(0.05, 1))

fig_1 <- plot_grid(
  diagram_1_plot,
  fig_dsp_img_1,
  NULL,
  fig_dsp_img,
  diagram_2_plot,
  ncol = 1,
  rel_heights = c(1, 1, 0.05, 0.8, 0.8),
  labels = c("A", "", "", "D", "E", "F")
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(
  "./Results/Figures/Final_figures/Figure_1.png",
  plot = fig_1,
  units = "cm",
  width = 21,
  height = 26,
  dpi = 600
)
