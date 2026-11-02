# Called after panel analyses; write figures only, never diagnostic tables.
source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
renderer <- make_renderer(analysis)
for (panel in panels) {
  plot_files <- c(plot_files, save_panel(renderer, panel, output_dir, output_format, dpi))
}
if (setequal(panels, LETTERS[1:5])) {
  plot_files <- c(plot_files, save_panel(renderer, "complete", output_dir, output_format, dpi))
}
cat("Fig2: PDF/PNG output in ", output_dir, "\n", sep = "")
