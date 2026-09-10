# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

root <- file.path(fig3_work_dir(), "Fig2/Newman_Monaco_7tools_validation")
worker <- fig3_script("pbmc/compare_external_pbmc_7tools_worker.R")
n_cores <- fig3_cores()

datasets <- list(
  list(
    dataset_name = "RNA_Sieve_newman_pbmcs",
    immucellai2_mode = "flat",
    immucellai2_add_unknown = FALSE
  ),
  list(
    dataset_name = "RNA_Sieve_monaco_pbmcs",
    immucellai2_mode = "tcell",
    immucellai2_add_unknown = TRUE
  )
)

requested <- if (exists("fig3_requested_datasets", inherits = TRUE)) {
  fig3_requested_datasets
} else commandArgs(trailingOnly = TRUE)
if (length(requested)) {
  datasets <- Filter(function(x) x$dataset_name %in% requested, datasets)
  if (!length(datasets)) stop("Requested dataset was not recognized.")
}

for (cfg in datasets) {
  dataset_dir <- file.path(root, cfg$dataset_name)
  std_dir <- file.path(dataset_dir, "standardized")
  out_dir <- file.path(dataset_dir, "seven_tools_compare")

  e <- new.env(parent = globalenv())
  e$dataset_name <- cfg$dataset_name
  e$bulk_file <- file.path(std_dir, "expression.csv")
  e$truth_file <- file.path(std_dir, "true_proportions.csv")
  e$out_dir <- out_dir
  e$n_cores <- n_cores
  e$chain_length <- 600
  e$burn_in <- 500
  e$immucellai2_mode <- cfg$immucellai2_mode
  e$immucellai2_add_unknown <- cfg$immucellai2_add_unknown
  e$immucellai2_method_name <- paste0(
    "ImmuCellAI2_", cfg$immucellai2_mode, "_VB",
    if (cfg$immucellai2_add_unknown) "_UNKNOWN" else ""
  )

  message("\n========== Running ", cfg$dataset_name, " ==========")
  sys.source(worker, envir = e)
  if (exists("run_status", envir = e, inherits = FALSE) &&
      any(e$run_status$Status == "FAILED")) stop("A method failed in ", cfg$dataset_name)
}

all_dataset_names <- c("RNA_Sieve_newman_pbmcs", "RNA_Sieve_monaco_pbmcs")
summary_files <- file.path(root, all_dataset_names, "seven_tools_compare",
                           "comparison_summary_7tools.txt")
available <- file.exists(summary_files)
combined <- do.call(rbind, lapply(which(available), function(i) {
  x <- read.delim(summary_files[i], check.names = FALSE)
  x$Dataset <- all_dataset_names[i]
  x
}))
combined <- combined[, c("Dataset", setdiff(names(combined), "Dataset"))]
write.table(
  combined,
  file.path(root, "Newman_Monaco_7tools_comparison_summary.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

best <- do.call(rbind, lapply(split(combined, combined$Dataset), function(x) {
  x[which.max(x$PenalizedMeanPearson), , drop = FALSE]
}))
rownames(best) <- NULL
write.table(
  best,
  file.path(root, "Newman_Monaco_7tools_best_method.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(best)
