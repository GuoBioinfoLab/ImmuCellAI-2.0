require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required R packages missing: ", paste(missing, collapse = ", "))
  for (package in packages) suppressPackageStartupMessages(library(package, character.only = TRUE))
}
require_packages(c("ggplot2", "dplyr", "cowplot"))

# Preserve the original density, boxplot and significance-label definitions.
sig_symbol <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

make_half_violin_tb <- function(dat, width = 0.25) {
  dat <- dat %>% filter(!is.na(Abundance), !is.na(status), !is.na(CellType))
  status_levels <- levels(dat$status)
  pieces <- list()
  k <- 1L
  for (status_i in status_levels) {
    vals <- dat$Abundance[dat$status == status_i]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) next
    dens <- density(vals, n = 128, from = min(vals), to = max(vals), na.rm = TRUE)
    dens_scaled <- dens$y / max(dens$y, na.rm = TRUE) * width
    x_base <- ifelse(status_i == "Active TB", 1.20, 2.20)
    pieces[[k]] <- data.frame(
      status = factor(status_i, levels = status_levels),
      x = c(x_base + dens_scaled, rep(x_base, length(dens$x))),
      y = c(dens$x, rev(dens$x)),
      group_id = status_i,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  bind_rows(pieces)
}

make_tb_panel <- function(cell, dat, stat_df, show_y = FALSE) {
  tmp <- dat %>%
    filter(CellType == cell) %>%
    mutate(
      x = ifelse(status == "Active TB", 1, 2),
      x_box = x - 0.16
    )
  half_df <- make_half_violin_tb(tmp)
  pval <- stat_df$P_Wilcox[as.character(stat_df$CellType) == cell]
  y_max <- max(tmp$Abundance, na.rm = TRUE)
  y_min <- min(tmp$Abundance, na.rm = TRUE)
  y_pad <- ifelse(y_max > y_min, (y_max - y_min) * 0.18, max(y_max, 1e-6) * 0.25)
  y_label <- y_max + y_pad * 0.45

  ggplot(tmp, aes(x = x_box, y = Abundance)) +
    geom_polygon(
      data = half_df,
      aes(x = x, y = y, group = group_id, fill = status),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8,
      alpha = 0.9
    ) +
    stat_boxplot(
      aes(group = status, colour = status),
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.75
    ) +
    geom_boxplot(
      aes(group = status, colour = status, fill = status),
      width = 0.34,
      outlier.shape = NA,
      linewidth = 0.75,
      alpha = 0.95
    ) +
    stat_summary(
      aes(group = status),
      fun = "median",
      geom = "point",
      shape = 16,
      size = 2.1,
      color = "white"
    ) +
    annotate("text", x = 1.5, y = y_label, label = sig_symbol(pval),
             size = 8, fontface = "bold", color = "black") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Active TB", "Latent TB"),
      limits = c(0.48, 2.62),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c("Active TB" = "#E15759FF", "Latent TB" = "#4E79A7FF")) +
    scale_color_manual(values = c("Active TB" = "black", "Latent TB" = "black")) +
    labs(x = NULL, y = if (show_y) "Cell Ratio" else NULL, title = cell) +
    coord_cartesian(ylim = c(y_min, y_max + y_pad), clip = "off") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", color = "#8B1A1A", hjust = 0),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(size = 12, color = "#8B1A1A", face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.8, linetype = "dashed"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.05),
      aspect.ratio = 1,
      plot.margin = margin(4, 8, 4, 8)
    )
}

make_half_violin_data <- function(dat, x_map, width = 0.25) {
  dat <- dat %>% filter(!is.na(Abundance), !is.na(status), !is.na(CellType))
  pieces <- list()
  k <- 1L
  for (status_i in names(x_map)) {
    vals <- dat$Abundance[as.character(dat$status) == status_i]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) next
    dens <- density(vals, n = 128, from = min(vals), to = max(vals), na.rm = TRUE)
    dens_scaled <- dens$y / max(dens$y, na.rm = TRUE) * width
    x_base <- unname(x_map[status_i]) + 0.20
    pieces[[k]] <- data.frame(
      status = factor(status_i, levels = names(x_map)),
      x = c(x_base + dens_scaled, rep(x_base, length(dens$x))),
      y = c(dens$x, rev(dens$x)),
      group_id = status_i,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  bind_rows(pieces)
}

make_cell_panel <- function(dat, stat_df, cell, comparisons, x_map, show_y = FALSE) {
  tmp <- dat %>%
    filter(CellType == cell) %>%
    mutate(
      x = unname(x_map[as.character(status)]),
      x_box = x - 0.16
    )
  half_df <- make_half_violin_data(tmp, x_map = x_map)
  y_max <- max(tmp$Abundance, na.rm = TRUE)
  y_min <- min(tmp$Abundance, na.rm = TRUE)
  y_range <- ifelse(y_max > y_min, y_max - y_min, max(y_max, 1e-6))
  y_pad <- y_range * 0.22
  stat_plot <- stat_df %>%
    filter(CellType == cell) %>%
    mutate(
      x = unname(x_map[group1]),
      xend = unname(x_map[group2]),
      y = y_max + seq_along(p) * y_range * 0.10,
      yend = y,
      label_y = y + y_range * 0.035
    )

  ggplot(tmp, aes(x = x_box, y = Abundance)) +
    geom_polygon(
      data = half_df,
      aes(x = x, y = y, group = group_id, fill = status),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.75,
      alpha = 0.9
    ) +
    stat_boxplot(
      aes(group = status, colour = status),
      geom = "errorbar",
      width = 0.16,
      linewidth = 0.7
    ) +
    geom_boxplot(
      aes(group = status, colour = status, fill = status),
      width = 0.32,
      outlier.shape = NA,
      linewidth = 0.7,
      alpha = 0.95
    ) +
    stat_summary(
      aes(group = status),
      fun = "median",
      geom = "point",
      shape = 16,
      size = 1.8,
      color = "white"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.65
    ) +
    geom_text(
      data = stat_plot,
      aes(x = (x + xend) / 2, y = label_y, label = p.signif),
      inherit.aes = FALSE,
      size = 5.5,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = unname(x_map),
      labels = names(x_map),
      limits = c(min(unname(x_map)) - 0.55, max(unname(x_map)) + 0.65),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c(
      "Died" = "#F28E2BFF",
      "Survived" = "#76B7B2FF",
      "Died-IRIS" = "#F28E2BFF",
      "Survived-IRIS" = "#76B7B2FF",
      "Early ART" = "#00A087FF",
      "Deferred ART" = "#EFC000FF"
    )) +
    scale_color_manual(values = setNames(rep("black", length(x_map)), names(x_map))) +
    labs(x = NULL, y = if (show_y) "Cell Ratio" else NULL, title = cell) +
    coord_cartesian(ylim = c(y_min, y_max + y_pad + nrow(stat_plot) * y_range * 0.08), clip = "off") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", color = "#8B1A1A", hjust = 0),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(size = 12, color = "#8B1A1A", face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75, linetype = "dashed"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.95),
      aspect.ratio = 1,
      plot.margin = margin(4, 8, 12, 8)
    )
}


prepare_plot_data <- function(panel) {
  result <- analysis[[panel]]
  result$plot_data$CellType <- factor(result$plot_data$CellType, levels = result$cells)
  result$plot_data$status <- factor(result$plot_data$status, levels = result$groups)
  result
}

save_plot <- function(plot, panel) {
  files <- character()
  for (extension in if (output_format == "both") c("pdf", "png") else output_format) {
    path <- file.path(output_dir, paste0("Fig6_", panel, ".", extension))
    ggplot2::ggsave(path, plot, width = 12.5, height = 3.2, dpi = dpi, bg = "white")
    files <- c(files, path)
  }
  files
}

combine_box_plots <- function(plots, panel) {
  cowplot::plot_grid(plotlist = plots, nrow = 1, align = "hv",
                    labels = c(panel, "", "", ""), label_size = 18,
                    label_fontface = "bold", hjust = -0.2, vjust = 1.4)
}

plot_comparison_panel <- function(panel) {
  input <- prepare_plot_data(panel)
  x_map <- setNames(seq_along(input$groups), input$groups)
  plots <- lapply(seq_along(input$cells), function(i) {
    plot <- make_cell_panel(input$plot_data, input$statistics, input$cells[i],
                            list(input$groups), x_map, i == 1)
    # Match the existing archive's composite-label layout.
    plot$layers <- Filter(function(layer) !inherits(layer$geom, "GeomSegment"), plot$layers)
    plot + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  })
  save_plot(combine_box_plots(plots, panel), panel)
}
