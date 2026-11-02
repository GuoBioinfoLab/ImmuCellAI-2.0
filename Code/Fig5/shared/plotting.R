# Plot definitions extracted from the original Fig5 scripts without changing styles.
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required R packages missing: ", paste(missing, collapse = ", "))
  for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
}
require_packages(c("ggplot2", "dplyr", "cowplot", "ComplexHeatmap", "grid"))
theme_blue <- theme(
  plot.title = element_text(size = 13, face = "bold", color = "darkred", hjust = 0, lineheight = 1.2),
  plot.subtitle = element_text(size = 11, face = "bold", color = "grey30", lineheight = 1.2, hjust = 0),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "gray80", linewidth = 0.7, linetype = "dashed"),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(vjust = 1, face = "bold", size = 12, color = "darkred"),
  axis.title.y = element_text(size = 12, face = "bold", color = "darkred"),
  axis.text.x = element_text(size = 10, colour = "black"),
  axis.text.y = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 11, colour = "black"),
  legend.text = element_text(size = 10, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
  legend.key = element_blank(),
  strip.background = element_rect(fill = "#F7E6E8", color = "black", linewidth = 0.6),
  strip.text = element_text(size = 10, colour = "black", face = "bold")
)
age_cols <- c(
  "age0-1" = "#E15759FF",
  "age10-20" = "#4E79A7FF",
  "age20-30" = "#F28E2BFF",
  "age30-50" = "#76B7B2FF",
  "age50-70" = "#59A14FFF",
  "age70+" = "#EDC948FF"
)
cell_y_limits <- list(
  M0 = c(0, 0.018),
  Bn = c(0, 0.01),
  CD8Tem = c(0, 0.0038)
)
make_cell_plot <- function(cell, plot_data, jitter_size = 1.8) {
  dat <- plot_data %>% filter(cell_type == cell)
  p <- ggplot(dat, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.8) +
    geom_jitter(shape = 21, color = "gray2", position = position_jitter(0.17), size = jitter_size, alpha = 0.78) +
    scale_x_discrete(labels = c("0-1", "10-20", "20-30", "30-50", "50-70", "70+")) +
    scale_fill_manual(values = age_cols, drop = FALSE) +
    labs(title = cell, x = "Age Group", y = "Cell Ratio", fill = NULL) +
    theme_blue +
    geom_smooth(
      method = "loess", formula = y ~ x, linewidth = 1, se = FALSE,
      linetype = "solid", aes(group = 1), color = "#3F88C5"
    ) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  if (cell %in% names(cell_y_limits)) {
    p <- p + coord_cartesian(ylim = cell_y_limits[[cell]])
  }
  p
}
save_plot <- function(plot, panel, width, height, heatmap = FALSE) {
  files <- character()
  draw_heatmap <- function(path, ext) {
    if (ext == "pdf") cairo_pdf(path, width = width, height = height)
    else png(path, width = width, height = height, units = "in", res = dpi, type = "cairo")
    on.exit(dev.off(), add = TRUE)
    ComplexHeatmap::draw(plot, heatmap_legend_side = "left", annotation_legend_side = "left",
                         newpage = TRUE, merge_legend = TRUE)
  }
  for (ext in if (output_format == "both") c("pdf", "png") else output_format) {
    path <- file.path(output_dir, paste0("Fig5_", panel, ".", ext))
    set.seed(123)
    if (heatmap) draw_heatmap(path, ext)
    else ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
    files <- c(files, path)
  }
  files
}
plot_box_panel <- function(panel) {
  plot_data <- analysis[[panel]]
  plot_data$group <- factor(plot_data$group, levels = age_levels)
  plot <- cowplot::plot_grid(
    plotlist = lapply(panel_cells[[panel]], make_cell_plot, plot_data = plot_data),
    ncol = if (panel == "C") 3 else 4, align = "hv"
  )
  save_plot(plot, panel, width = if (panel == "C") 12 else 16,
            height = if (panel == "C") 13.2 else 4.4)
}
