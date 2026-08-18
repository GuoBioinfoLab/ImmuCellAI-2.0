args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
scripts <- c(
  "01_deconvolve_tcga.R",
  "02_consensus_clustering.R",
  "03_plot_phenotypes_and_survival.R",
  "05_associate_score_with_tcga_phenotypes.R"
)
base <- "analysis/tumor_immune_phenotypes_immuicbscore"
for (script in scripts) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(file.path(base, script), config_file))
  if (status != 0L) stop(script, " failed with status ", status)
}
source(config_file)
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(base, "04_immuicbscore_pooled_fivefold_cv.R"),
    icb_fraction_file, icb_clinical_file, icb_cv_output_dir, as.character(n_cores))
)
if (status != 0L) stop("ImmuICBscore cross-validation failed.")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(base, "06_plot_icb_response_cell_fractions.R"), config_file)
)
if (status != 0L) stop("ICB response cell-fraction plotting failed.")
