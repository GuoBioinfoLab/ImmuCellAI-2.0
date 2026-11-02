source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
plots <- list()
for (panel in panels) {
  plots[[panel]] <- build_panel(panel)
  p <- plots[[panel]]
  files <- c(files, save_plot(p$plot, paste0("Supplement_Fig6_1022_", panel), p$width, p$height))
}
if (all(LETTERS[1:10] %in% panels)) files <- c(files, save_plot(combine_panels(plots[LETTERS[1:10]]), "Supplement_Fig6_1022_complete", 17, 23))
