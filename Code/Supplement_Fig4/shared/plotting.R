required <- c("ggplot2", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
save_plot <- function(plot, stem, width, height) {
  result <- character()
  for (extension in if (output_format == "both") c("pdf", "png") else output_format) {
    path <- file.path(output_dir, paste0(stem, ".", extension))
    if (extension == "pdf") {
      grDevices::cairo_pdf(path, width = width, height = height, family = "sans", bg = "white")
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = dpi,
                     type = "cairo", family = "sans", bg = "white")
    }
    tryCatch(print(plot), finally = invisible(grDevices::dev.off()))
    result <- c(result, path)
  }
  result
}

build_panel <- function(panel) {
  result <- analysis[[panel]]; meta <- result$meta; curve <- result$curve
  regimen <- if (isTRUE(getOption("Supplement_Fig4.verified_labels", FALSE))) meta$verified_regimen else meta$original_figure_regimen
  curve$heading <- paste(meta$figure_cancer_label, paste(meta$accession, regimen),
                         sprintf("AUC = %.3f, N = %d", result$AUC, meta$N), sep = "\n")
  p <- ggplot(curve, aes(FPR, TPR)) +
    geom_abline(slope = 1, intercept = 0, color = "#909090", linetype = "dashed", linewidth = .7) +
    geom_path(color = "#379DA5", linewidth = 1.1) +
    scale_x_continuous(breaks = c(.25, .5, .75, 1), labels = function(x) sprintf("%.2f", x), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, .25), labels = function(x) sprintf("%.2f", x), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    facet_wrap(~heading) +
    labs(x = "False positive rate", y = "True positive rate") +
    theme_bw(base_size = 11, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#DDDDDD", linewidth = .45),
          panel.border = element_rect(color = "black", linewidth = .7),
          strip.background = element_rect(fill = "#F7E7EA", color = "black", linewidth = .7),
          strip.text = element_text(face = "bold", size = 10.5, margin = margin(5, 3, 5, 3)),
          axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
          aspect.ratio = 1, plot.margin = margin(5, 15, 5, 5))
  list(plot = p, width = 4.6, height = 5.3)
}
combine_panels <- function(plots) {
  display <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]]$plot + labs(x = NULL, y = NULL)
    if (i <= 4) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    if (!i %in% c(1, 5)) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p
  })
  body <- plot_grid(plotlist = display, ncol = 4, align = "hv")
  ggdraw() + draw_plot(body, x = .04, y = .05, width = .955, height = .90) +
    draw_label("Per-cohort ROC curves", x = .055, y = .99, hjust = 0, vjust = 1, fontface = "bold", size = 16) +
    draw_label("False positive rate", x = .52, y = .012, fontface = "bold", size = 12) +
    draw_label("True positive rate", x = .012, y = .50, angle = 90, fontface = "bold", size = 12)
}
