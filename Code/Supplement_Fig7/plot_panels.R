source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
plots <- list()
for (panel in panels) {
  plots[[panel]] <- build_panel(panel)
  result <- plots[[panel]]
  files <- c(files, save_plot(result$plot, paste0("Supplement_Fig7_", panel), result$width, result$height))
}
if (all(c("A", "B") %in% panels)) {
  combined <- plot_grid(plots$A$body, plots$B$body, plots$B$legend, ncol = 1, rel_heights = c(1, 2, .12))
  files <- c(files, save_plot(combined, "Supplement_Fig7_complete", 14, 11.4))
}
