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

scenario_order <- c("B", "DC", "Monocyte", "TNK")
scenario_titles <- c("B-cell reference omitted", "DC reference omitted", "Monocyte reference omitted", "T/NK-cell reference omitted")
base_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(panel.grid.major = element_line(color = "#E4E4E4", linewidth = .35),
        panel.grid.minor = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = .5, size = 12),
        plot.margin = margin(10, 12, 8, 8), legend.position = "none")
build_panel <- function(panel) {
  ps <- lapply(seq_along(scenario_order), function(i) {
    scenario <- scenario_order[i]
    if (panel == "A") {
      d <- analysis$A$points; d <- d[d$scenario == scenario, ]
      m <- analysis$A$metrics; m <- m[m$scenario == scenario, ]
      p <- ggplot(d, aes(omitted_truth_percent, UNKNOWN_weight_percent)) +
        geom_point(shape = 21, fill = "#376EA8", color = "white", stroke = .3, size = 2.6, alpha = .9) +
        geom_segment(x = min(d$omitted_truth_percent), xend = max(d$omitted_truth_percent),
          y = 100 * m$intercept + m$slope * min(d$omitted_truth_percent),
          yend = 100 * m$intercept + m$slope * max(d$omitted_truth_percent),
          inherit.aes = FALSE, color = "#376EA8", linewidth = .8) +
        annotate("text", x = max(d$omitted_truth_percent) * .054, y = max(d$UNKNOWN_weight_percent) * 1.35,
                 label = sprintf("Pearson r = %.3f\nSlope = %.3f", m$r, m$slope),
                 hjust = 0, vjust = 1, size = 3.8) +
        scale_x_continuous(limits = c(0, max(d$omitted_truth_percent) * 1.08), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, max(d$UNKNOWN_weight_percent) * 1.38), expand = c(0, 0)) +
        labs(x = "Omitted-cell fraction by flow (%)", y = if (i == 1) "Estimated UNKNOWN weight (%)" else NULL)
    } else {
      d <- analysis$B[analysis$B$scenario == scenario, ]; d <- d[match(c("Incomplete / OFF", "Incomplete / ON"), d$condition), ]
      d$state <- factor(c("OFF", "ON"), c("OFF", "ON"))
      p <- ggplot(d, aes(state, mae_pp, fill = state)) +
        geom_col(width = .52, color = "#686868", linewidth = .5) +
        geom_errorbar(aes(ymin = ci95_low_pp, ymax = ci95_high_pp), width = .16, color = "#FEC466", linewidth = .8) +
        geom_text(aes(y = ci95_high_pp + .55, label = sprintf("%.2f", mae_pp)), vjust = 0, size = 3.8) +
        scale_fill_manual(values = c(OFF = "#A2A2A2", ON = "#376EA8")) +
        scale_y_continuous(limits = c(0, max(d$ci95_high_pp) * 1.24 + .5), expand = c(0, 0)) +
        labs(x = NULL, y = if (i == 1) "MAE (percentage points)" else NULL)
    }
    p + labs(title = scenario_titles[i]) + base_theme
  })
  heading <- if (panel == "A") "UNKNOWN weight versus omitted-cell fraction" else "Estimation error for remaining known immune groups"
  body <- plot_grid(plotlist = ps, nrow = 1, align = "h")
  p <- ggdraw() + draw_plot(body, y = .10, height = .72) +
    draw_label(panel, x = .006, y = .95, hjust = 0, fontface = "bold", size = 21) +
    draw_label(heading, x = .06, y = .95, hjust = 0, fontface = "bold", size = 15)
  if (panel == "B") p <- p + draw_label("UNKNOWN OFF", x = .39, y = .035, hjust = 0, size = 12) +
    draw_label("UNKNOWN ON", x = .55, y = .035, hjust = 0, size = 12) +
    draw_grob(grid::rectGrob(gp = grid::gpar(fill = "#A2A2A2", col = "#686868")), x = .367, y = .021, width = .014, height = .028) +
    draw_grob(grid::rectGrob(gp = grid::gpar(fill = "#376EA8", col = "#686868")), x = .527, y = .021, width = .014, height = .028)
  list(plot = p, body = body, width = 16, height = 4.8)
}
combine_panels <- function(plots) ggdraw() +
  draw_plot(plots$A$plot, y = .48, height = .48) + draw_plot(plots$B$plot, y = 0, height = .48) +
  draw_label("Real-bulk reference-omission validation", x = .025, y = .995, hjust = 0, vjust = 1, fontface = "bold", size = 19)
