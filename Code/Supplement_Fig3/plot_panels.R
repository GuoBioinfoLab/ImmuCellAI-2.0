source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
plots <- list()
for (panel in panels) {
  plots[[panel]] <- build_panel(panel)
  p <- plots[[panel]]
  files <- c(files, save_plot(p$plot, paste0("Supplement_Fig3_", panel), p$width, p$height))
}
if (all(c("A", "B") %in% panels)) files <- c(files, save_plot(combine_panels(plots[c("A", "B")]), "Supplement_Fig3_complete", 16, 8.6))
