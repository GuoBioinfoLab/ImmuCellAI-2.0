args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/infectious_disease_states/config.R"
base <- "analysis/infectious_disease_states"
for (script in c("01_prepare_and_deconvolve.R", "02_plot_active_vs_latent_tb.R", "03_plot_hiv_clinical_comparisons.R", "04_hiv_longitudinal_interaction.R")) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(file.path(base, script), config_file))
  if (status != 0L) stop(script, " failed with status ", status)
}
