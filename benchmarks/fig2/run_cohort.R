args <- commandArgs(trailingOnly = TRUE)
cohorts <- c("GSE164522", "GSE176078", "GSE146771")
if (!length(args) || !args[1] %in% c(cohorts, "all", "plot")) {
  stop("Usage: Rscript benchmarks/fig2/run_cohort.R ",
       "{GSE164522|GSE176078|GSE146771|all|plot} [config.R]")
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
script_dir <- dirname(script_file)
repo_root <- dirname(dirname(script_dir))
config_file <- if (length(args) >= 2L) normalizePath(args[2], mustWork = TRUE) else NULL
setwd(repo_root)

defaults <- new.env(parent = baseenv())
sys.source(file.path(script_dir, "config.example.R"), envir = defaults)
config <- defaults$config
if (!is.null(config_file)) {
  user_config <- new.env(parent = baseenv())
  sys.source(config_file, envir = user_config)
  if (!exists("config", envir = user_config, inherits = FALSE) || !is.list(user_config$config)) {
    stop("The configuration file must define a list named config.")
  }
  extra <- setdiff(names(user_config$config), names(config))
  if (length(extra)) stop("Unknown configuration field(s): ", paste(extra, collapse = ", "))
  config <- utils::modifyList(config, user_config$config, keep.null = TRUE)
}
if (length(config$library_paths)) {
  if (any(!dir.exists(config$library_paths))) stop("An R library path does not exist.")
  .libPaths(c(config$library_paths, .libPaths()))
}
if (length(config$n_cores) != 1L || !is.finite(config$n_cores) || config$n_cores < 1) {
  stop("n_cores must be a positive integer.")
}
if (config$burn_in < 0 || config$chain_length <= config$burn_in) {
  stop("chain_length must exceed burn_in, and burn_in must be non-negative.")
}

run_one <- function(dataset) {
  env <- list2env(config[setdiff(names(config), c("output_dir", "bulk_files", "library_paths"))],
                  parent = globalenv())
  # NULL leaves the original script's packaged-resource default in effect.
  for (key in c("reference_file", "immune_gene_file")) {
    if (is.null(env[[key]])) rm(list = key, envir = env)
  }
  if (!is.null(config$bulk_files) && dataset %in% names(config$bulk_files)) {
    env$bulk_file <- unname(config$bulk_files[[dataset]])
  }
  env$out_dir <- file.path(config$output_dir, dataset)
  dir.create(env$out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(script_dir, paste0("compare_", tolower(dataset),
                                     "_variable_background_7tools.R"))
  on.exit(writeLines(capture.output(sessionInfo()), file.path(env$out_dir, "sessionInfo.txt")))
  saveRDS(config, file.path(env$out_dir, "run_config.rds"))
  sys.source(path, envir = env)
  inputs <- c(script = path, bulk = env$bulk_file, reference = env$reference_file,
              markers = env$immune_gene_file)
  write.table(data.frame(Role = names(inputs), File = unname(inputs),
                         MD5 = unname(tools::md5sum(inputs))),
              file.path(env$out_dir, "input_checksums.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  if (any(env$run_status$Status == "FAILED")) {
    stop(dataset, ": one or more methods failed; inspect run_status.txt before plotting.")
  }
  invisible(NULL)
}

if (args[1] == "all") {
  for (dataset in cohorts) run_one(dataset)
} else if (args[1] != "plot") {
  run_one(args[1])
}
if (args[1] %in% c("all", "plot")) {
  env <- new.env(parent = globalenv())
  env$out_dir <- file.path(config$output_dir, "plots")
  env$input_files <- data.frame(
    Dataset = cohorts,
    File = file.path(config$output_dir, cohorts, "all_methods_per_celltype_metrics.txt")
  )
  sys.source(file.path(script_dir,
                      "plot_4datasets_7tools_correlation_pie_matrix_final_ordered.R"), envir = env)
}
