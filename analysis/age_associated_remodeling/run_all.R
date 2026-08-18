args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
base <- "analysis/age_associated_remodeling"
for (script in c("01_prepare_and_deconvolve.R", "02_plot_age_cell_abundance.R", "03_plot_age_pattern_heatmap.R", "04_plot_major_lineage_composition.R")) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(file.path(base, script), config_file))
  if (status != 0L) stop(script, " failed with status ", status)
}
