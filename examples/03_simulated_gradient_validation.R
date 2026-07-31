args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript examples/03_simulated_gradient_validation.R",
      "prediction_file truth_file output_dir [mapping_file]"
    )
  )
}

prediction_file <- normalizePath(args[1], mustWork = TRUE)
truth_file <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
mapping_file <- if (length(args) >= 4L) normalizePath(args[4], mustWork = TRUE) else NA_character_

script <- file.path("benchmarks", "evaluate_predictions.R")
if (!file.exists(script)) {
  stop("Run this command from the repository root; missing ", script)
}

command <- c(script, prediction_file, truth_file, output_dir)
if (!is.na(mapping_file)) command <- c(command, mapping_file)

status <- system2(file.path(R.home("bin"), "Rscript"), command)
if (!identical(status, 0L)) {
  stop("Prediction evaluation failed with exit status ", status)
}
