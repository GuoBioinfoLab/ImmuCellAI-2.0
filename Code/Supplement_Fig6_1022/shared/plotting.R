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

titles <- c(A = "Validation-set permutation importance", B = "Top-10 inclusion frequency",
 C = "ImmuICBscore", D = "Major immune groups (RF)", E = "CD45 expression",
 F = "ImmuICBscore", G = "Major immune groups (RF)", H = "CD45 expression",
 I = "ICB-response discrimination", J = "Adjusted survival association of ImmuICBscore")
base_theme <- theme_classic(base_size = 10, base_family = "sans") +
 theme(axis.text = element_text(color = "black"), plot.margin = margin(12, 18, 12, 12))
decorate <- function(p, panel) ggdraw() + draw_plot(p, y = .01, height = .87) +
 draw_label(panel, x = .015, y = .965, hjust = 0, fontface = "bold", size = 19) +
 draw_label(titles[panel], x = .13, y = .965, hjust = 0, fontface = "bold", size = 12)
build_panel <- function(panel) {
 if (panel %in% c("A", "B")) {
   d <- analysis[[panel]]; d$feature <- factor(d$feature, levels = rev(d$feature))
   if (panel == "A") {
     p <- ggplot(d, aes(mean_auc_decrease, feature)) +
       geom_vline(xintercept = 0, color = "#AAAAAA", linetype = "dashed") +
       geom_errorbar(aes(xmin = mean_auc_decrease - sd_auc_decrease, xmax = mean_auc_decrease + sd_auc_decrease),
                     orientation = "y", width = .22, color = "#A78E92") +
       geom_point(color = "#9B3545", size = 2) +
       scale_x_continuous(breaks = seq(0, .0175, .0025)) +
       coord_cartesian(xlim = c(-.0008, .0177)) +
       labs(x = "Decrease in validation AUC (mean +/- SD)", y = NULL)
   } else {
     p <- ggplot(d, aes(100 * top10_frequency, feature)) + geom_col(fill = "#287E79", width = .62) +
       geom_text(aes(x = 100 * top10_frequency + 1.5, label = paste0(top10_count, "/20")), hjust = 0, size = 3.2) +
       scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25), expand = c(0, 0)) +
       labs(x = "Repeats with a top-10 rank (%)", y = NULL)
   }
   p <- p + base_theme + theme(panel.grid.major.x = element_line(color = "#EDEDED"), axis.ticks.y = element_blank())
   width <- 8; height <- 6.5
 } else if (panel %in% LETTERS[3:8]) {
   result <- analysis[[panel]]; d <- result$curve; row <- result$statistics
   starts <- d[match(c("High", "Low"), d$group), ]
   starts$months <- 0; starts$survival <- starts$lower <- starts$upper <- 1; starts$n_censor <- 0
   d <- rbind(starts, d); d <- d[order(d$group, d$months), ]
   # Duplicate each event time to draw stepwise confidence bands, matching the KM estimator.
   bands <- do.call(rbind, lapply(split(d, d$group), function(g) {
     z <- g[rep(seq_len(nrow(g)), each = 2), ]
     z$months <- c(0, as.vector(rbind(tail(g$months, -1), tail(g$months, -1))), tail(g$months, 1))
     z
   }))
   pooled <- panel %in% LETTERS[3:5]
   label <- if (pooled) "pooled" else "cancer-specific"
   note <- if (pooled) sprintf("Log-rank P = %.3g", row$logrank_p) else
     sprintf("Log-rank P = %.3g\nCancer-stratified P = %.3g\nAdjusted HR = %.3f (%.3f-%.3f)",
       row$logrank_p, row$cancer_stratified_logrank_p, row$adjusted_HR, row$adjusted_CI_low, row$adjusted_CI_high)
   p <- ggplot(d, aes(months, survival, color = group, group = group)) +
     geom_ribbon(data = bands[is.finite(bands$lower) & is.finite(bands$upper), ], aes(ymin = lower, ymax = upper, fill = group), alpha = .10, color = NA) +
     geom_step(linewidth = .6) +
     geom_point(data = d[d$n_censor > 0, ], shape = 3, size = .7, stroke = .25, alpha = .6) +
     scale_color_manual(values = c(High = "#A23649", Low = "#35817D"),
       labels = c(sprintf("High (n=%s)", format(row$high_n, big.mark = ",")), sprintf("Low (n=%s)", format(row$low_n, big.mark = ",")))) +
     scale_fill_manual(values = c(High = "#A23649", Low = "#35817D"), guide = "none") +
     scale_x_continuous(breaks = seq(0, 120, 24)) + scale_y_continuous(breaks = seq(0, 1, .25)) +
     coord_cartesian(xlim = c(0, 120), ylim = c(0, 1.02), expand = FALSE) +
     annotate("text", x = 4.2, y = .065, label = note, hjust = 0, vjust = 0, size = 3.0, lineheight = 1.3) +
     labs(x = "Months", y = if (panel %in% c("C", "F")) "Overall survival" else NULL, color = NULL,
          caption = paste0("High > ", label, " median; Low <= median")) +
     base_theme + theme(panel.grid.major.y = element_line(color = "#EEEEEE"),
       legend.position = "inside", legend.position.inside = c(.76, .9), legend.background = element_blank(),
       legend.key.height = grid::unit(.3, "cm"), legend.text = element_text(size = 9),
       plot.caption = element_text(hjust = .5, size = 9, margin = margin(t = 12)))
   width <- 6.8; height <- 5.7
 } else {
   is_i <- panel == "I"; d <- if (is_i) analysis$I else analysis$J$data; n <- nrow(d); d$y <- n:1
   d$value <- if (is_i) d$AUC else d$HR
   # Text columns and the forest axis share the same explicit row coordinates.
   xr <- if (is_i) c(.3, 1.01) else c(.82, 1.04)
   ticks <- if (is_i) c(.5, .75, 1) else c(.85, .9, .95, 1)
   xpos <- function(v) .44 + .23 * (v - xr[1]) / diff(xr)
   p <- ggdraw() +
     draw_line(x = c(.44, .67), y = c(.15, .15), color = "black", linewidth = .4) +
     draw_line(x = rep(xpos(if (is_i) .5 else 1), 2), y = c(.15, .87), color = "#AAAAAA", linetype = "dashed", linewidth = .5)
   for (tick in ticks) p <- p +
     draw_line(x = rep(xpos(tick), 2), y = c(.15, .14), color = "black", linewidth = .4) +
     draw_label(sprintf("%.2f", tick), x = xpos(tick), y = .11, size = 9)
   p <- p + draw_label(if (is_i) "ROC AUC" else "HR per 1-SD higher ImmuICBscore", x = .555, y = .04, size = 9)
   labels <- if (is_i) d$display_label else c("Clinical", "Clinical + CD45", "Clinical + major groups", "Clinical + major groups\n+ CD45")
   for (i in seq_len(n)) {
     yy <- .15 + .72 * (1 - (i - .5) / n)
     col <- if (is_i) "#287E79" else "#A23649"
     p <- p + draw_line(x = xpos(c(d$CI_low[i], d$CI_high[i])), y = c(yy, yy), color = col, linewidth = .7) +
       draw_line(x = rep(xpos(d$CI_low[i]), 2), y = yy + c(-.007, .007), color = col, linewidth = .7) +
       draw_line(x = rep(xpos(d$CI_high[i]), 2), y = yy + c(-.007, .007), color = col, linewidth = .7) +
       draw_grob(grid::pointsGrob(x = .5, y = .5, pch = 15, size = grid::unit(2, "mm"), gp = grid::gpar(col = col)),
         x = xpos(d$value[i]) - .006, y = yy - .006, width = .012, height = .012)
     p <- p + draw_label(labels[i], x = .015, y = yy, hjust = 0, size = 9)
     txt <- if (is_i) sprintf("%.3f\n[%.3f, %.3f]", d$value[i], d$CI_low[i], d$CI_high[i]) else
       sprintf("%.3f (%.3f, %.3f)\nP = %.3g", d$value[i], d$CI_low[i], d$CI_high[i], d$p[i])
     p <- p + draw_label(txt, x = .71, y = yy, hjust = 0, size = 9)
     if (is_i) p <- p + draw_label(d$n[i], x = .38, y = yy, size = 9)
   }
   p <- p + draw_label(if (is_i) "AUC [95% CI]" else "HR [95% CI]; Wald P", x = .71, y = .93, hjust = 0, fontface = "bold", size = 10)
   if (is_i) p <- p + draw_label("N", x = .38, y = .93, fontface = "bold", size = 10)
   width <- 11; height <- if (is_i) 7.2 else 5.7
 }
 list(plot = decorate(p, panel), width = width, height = height)
}
combine_panels <- function(plots) plot_grid(
 plot_grid(plots$A$plot, plots$B$plot, nrow = 1),
 plot_grid(plotlist = lapply(plots[LETTERS[3:5]], `[[`, "plot"), nrow = 1),
 plot_grid(plotlist = lapply(plots[LETTERS[6:8]], `[[`, "plot"), nrow = 1),
 plot_grid(plots$I$plot, plots$J$plot, nrow = 1), ncol = 1, rel_heights = c(1.06, 1, 1, 1.05))
