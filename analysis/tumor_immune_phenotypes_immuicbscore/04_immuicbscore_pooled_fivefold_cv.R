# This entry point keeps the case-study folder self-contained while using the
# single canonical implementation maintained in analysis/icb/.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript 04_immuicbscore_pooled_fivefold_cv.R fraction_file clinical_file output_dir [n_threads]")
}
canonical <- file.path("analysis", "icb", "immuicbscore_pooled_fivefold_cv.R")
status <- system2(file.path(R.home("bin"), "Rscript"), c(canonical, args))
if (status != 0L) stop("Canonical ImmuICBscore cross-validation script failed with status ", status)
