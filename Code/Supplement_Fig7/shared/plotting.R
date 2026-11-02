required <- c("ggplot2", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
save_plot <- function(plot, stem, width, height) {
  result <- character()
  for (extension in if (output_format == "both") c("pdf", "png") else output_format) {
    path <- file.path(output_dir, paste0(stem, ".", extension))
    if (extension == "pdf") {
      grDevices::cairo_pdf(path, width = width, height = height, family = "Arial", bg = "white")
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = dpi,
                     type = "cairo", family = "Arial", bg = "white")
    }
    tryCatch(print(plot), finally = invisible(grDevices::dev.off()))
    result <- c(result, path)
  }
  result
}

build_panel <- function(panel) {
  suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
  colors <- c(Died = "#A3C8DC", Survived = "#EA5D2D")
  burgundy <- "#8F2027"
  plots <- list()
  {
    cells <- read_tsv(file.path(data_dir, panel, "cell_order.tsv"))$CellType
    means <- read_tsv(file.path(data_dir, panel, "plot_data.tsv"))
    points <- read_tsv(file.path(data_dir, panel, "sample_level.tsv"))
    means$Time <- factor(means$Time, c("D0", "D1", "D4", "D8"))
    points$Time <- factor(points$Time, c("D0", "D1", "D4", "D8"))
    means$status <- factor(means$status, names(colors))
    points$status <- factor(points$status, names(colors))
    stats <- if (panel == "B") read_tsv(file.path(data_dir, panel, "statistics_displayed.tsv")) else NULL
    cell_plots <- lapply(seq_along(cells), function(i) {
      cell <- cells[i]
      summary <- means[means$CellType == cell, ]
      dat <- points[points$CellType == cell, ]
      p <- ggplot(summary, aes(Time, Mean, color = status, group = status))
      if (panel == "B") {
        p <- p + geom_point(data = dat, aes(y = Abundance),
                             position = position_jitter(width = .08, height = 0, seed = 7),
                             alpha = .18, size = 1.7, show.legend = FALSE) +
          geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = .12, linewidth = .6)
      }
      p <- p + geom_line(linewidth = .8) + geom_point(size = 2.8) +
        scale_color_manual(values = colors, drop = FALSE) +
        scale_x_discrete(drop = FALSE, expand = expansion(add = .5)) +
        labs(title = cell, x = "Time", y = if (i %% 4 == 1) "Cell Ratio" else NULL, color = NULL) +
        theme_bw(base_size = 11) +
        theme(strip.background = element_blank(),
              plot.title = element_text(size = 13, color = burgundy, face = "bold", hjust = 0),
              axis.text = element_text(color = "black", size = 11),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_text(color = burgundy, face = "bold", size = 12),
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              panel.grid.major.y = element_line(color = "grey80", linetype = "dashed", linewidth = .6),
              panel.border = element_rect(color = "black", fill = NA, linewidth = .85),
              legend.position = "bottom", legend.text = element_text(size = 11),
              plot.margin = margin(5, 10, 5, 5), aspect.ratio = 1)
      if (panel == "B") {
        row <- stats[stats$CellType == cell, ]
        limits <- range(c(dat$Abundance, summary$Mean - summary$SEM, summary$Mean + summary$SEM))
        span <- diff(limits)
        top <- limits[2] + .27 * span
        p <- p + annotate("text", x = 1, y = top,
                           label = sprintf("Pint (D0-D1) = %.3f\nFDR = %.3f", row$Pvalue, row$FDR),
                           hjust = 0, vjust = 1, size = 3.25, lineheight = 1.05) +
          scale_y_continuous(limits = c(limits[1] - .04 * span, top + .025 * span), expand = expansion(mult = 0))
      }
      p
    })
    legend <- get_legend(cell_plots[[1]])
    body <- plot_grid(plotlist = lapply(cell_plots, function(p) p + theme(legend.position = "none")), ncol = 4, align = "hv")
    labeled <- plot_grid(body, labels = panel, label_size = 17, label_fontface = "bold", label_x = 0, hjust = 0)
    plots[[panel]] <- labeled
    standalone <- plot_grid(labeled, legend, ncol = 1, rel_heights = c(if (panel == "A") 1 else 2, .12))
    list(plot = standalone, body = labeled, legend = legend, width = 14,
         height = if (panel == "A") 4.1 else 7.6)
  }
}
