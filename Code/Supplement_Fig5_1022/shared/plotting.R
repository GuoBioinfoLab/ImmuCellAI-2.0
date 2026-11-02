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

titles <- c(A = "Consensus cumulative distribution", B = "Relative change in CDF area",
            C = "Proportion of ambiguous clustering", D = "Immune phenotypes by cancer", E = "IP4 versus IP3: overall survival")
colors <- c("2" = "#616D79", "3" = "#B58B3A", "4" = "#A3374B", "5" = "#287D79", "6" = "#547BA6")
base_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = .35), axis.text = element_text(color = "black"),
        plot.margin = margin(10, 15, 10, 10), legend.background = element_blank())
decorate <- function(p, panel) ggdraw() + draw_plot(p, y = .02, height = .86) +
  draw_label(panel, x = .01, y = .96, hjust = 0, fontface = "bold", size = 20) +
  draw_label(titles[panel], x = .095, y = .96, hjust = 0, fontface = "bold", size = 12)
build_panel <- function(panel) {
  if (panel == "A") {
    d <- analysis$A$data
    p <- ggplot(d, aes(bin_midpoint, CDF, color = factor(K), group = K)) + geom_line(linewidth = .8) +
      scale_color_manual(values = colors, labels = paste("K =", 2:6), name = NULL) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.02), expand = FALSE) +
      labs(x = "Consensus index", y = "Cumulative probability") + base_theme +
      theme(legend.position = "inside", legend.position.inside = c(.85, .22), legend.key.height = grid::unit(.3, "cm"))
  } else if (panel %in% c("B", "C")) {
    d <- analysis[[panel]]; key <- if (panel == "B") "relative_delta_CDF_area" else "PAC_0.1_0.9"
    d$value <- d[[key]]; d <- d[!is.na(d$value), ]; d$color <- ifelse(d$K == 4, "#A3374B", "#6B7175")
    p <- ggplot(d, aes(K, value)) + geom_vline(xintercept = 4, color = "#A3374B", linetype = "dashed") +
      geom_line(color = "#6B7175", linewidth = .7) + geom_point(aes(color = color), size = 2.8) +
      geom_text(aes(label = sprintf("%.3f", value), color = color), vjust = -1.2, size = 3.6) +
      scale_color_identity() + scale_x_continuous(breaks = d$K) +
      coord_cartesian(xlim = if (panel == "B") c(2.7, 6.3) else c(1.7, 6.3), ylim = c(0, if (panel == "B") .48 else .36), expand = FALSE) +
      labs(x = "Number of clusters (K)", y = if (panel == "B") "Relative increase in CDF area" else "PAC (0.1 < consensus < 0.9)") + base_theme
  } else if (panel == "D") {
    d <- analysis$D$data
    p <- ggplot(d, aes(cancer, percent, fill = phenotype)) + geom_col(width = .75, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = c(IP1 = "#D59D46", IP2 = "#638EAF", IP3 = "#AD5266", IP4 = "#28877F"), name = NULL) +
      scale_x_discrete(labels = function(x) sub("^TCGA-", "", x)) +
      scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, .05))) +
      labs(x = NULL, y = "Patients (%)", caption = sprintf("N = 4,746; Cramer's V = %.3f; P < 0.001", analysis$D$cramers_V)) +
      base_theme + theme(legend.position = "top", axis.text.x = element_text(angle = 55, hjust = 1),
                        plot.caption = element_text(hjust = 0, size = 11, margin = margin(t = 15)))
  } else {
    d <- analysis$E$data; d$y <- 3:1
    p <- ggplot(d, aes(HR, y)) + geom_vline(xintercept = 1, color = "#999999", linetype = "dashed") +
      geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .13, color = "#237F7A", linewidth = .7) +
      geom_point(color = "#237F7A", shape = 15, size = 2.5) +
      scale_x_log10(limits = c(.35, 1.3), breaks = c(.4, .7, 1), labels = c("0.4", "0.7", "1.0")) +
      scale_y_continuous(limits = c(.5, 3.5), breaks = 3:1, labels = c("Unadjusted", "Cancer-stratified", "Cancer + clinical")) +
      labs(x = "HR (IP4 / IP3)", y = NULL) + base_theme + theme(panel.grid.major.y = element_blank(), panel.grid = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
    values <- ggplot(d, aes(0, y, label = sprintf("%.3f [%.3f, %.3f]\nP = %.4g", HR, CI_low, CI_high, p))) +
      geom_text(hjust = 0, size = 3.7) + scale_y_continuous(limits = c(.5, 3.5)) +
      coord_cartesian(xlim = c(0, 1), clip = "off") + theme_void() + labs(title = "HR [95% CI]") +
      theme(plot.title = element_text(size = 11, face = "bold"), plot.margin = margin(0, 0, 34, 0))
    body <- plot_grid(p, values, nrow = 1, rel_widths = c(.68, .32))
    pp <- analysis$E$comparisons$p[analysis$E$comparisons$model == "Clinical_adjusted"]
    p <- ggdraw() + draw_plot(body, y = .20, height = .78) +
      draw_label("Clinical adjustment: age, sex and stage.", x = .04, y = .12, hjust = 0, size = 11) +
      draw_label(sprintf("Overall adjusted IP1-IP4 test: P = %.5f.", pp), x = .04, y = .035, hjust = 0, size = 11)
  }
  list(plot = decorate(p, panel), width = if (panel == "D") 9 else if (panel == "E") 10.8 else 6.8,
       height = if (panel == "D") 5.7 else if (panel == "E") 5.4 else 4.8)
}
combine_panels <- function(plots) plot_grid(
  plot_grid(plotlist = lapply(plots[LETTERS[1:3]], `[[`, "plot"), nrow = 1),
  plot_grid(plots$D$plot, plots$E$plot, nrow = 1), ncol = 1, rel_heights = c(1, 1.05))
