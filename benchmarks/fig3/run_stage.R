args <- commandArgs(trailingOnly = TRUE)
stages <- c("prepare_pbmc", "compare_pbmc", "compare_xu", "xu_pairwise",
            "evaluate_xu", "plot_xu", "plot_pbmc", "runtime", "runtime_summary",
            "plot_runtime", "summarize", "compose", "plots", "all")
if (!length(args) || !args[1] %in% stages) {
  stop("Usage: Rscript benchmarks/fig3/run_stage.R {", paste(stages, collapse = "|"), "} [config.R]")
}
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), mustWork = TRUE)
home <- dirname(script_file)
repo <- dirname(dirname(home))
user_file <- if (length(args) > 1L) normalizePath(args[2], mustWork = TRUE) else NULL
setwd(repo)
defaults <- new.env(parent = baseenv())
sys.source(file.path(home, "config.example.R"), envir = defaults)
config <- defaults$config
if (!is.null(user_file)) {
  user <- new.env(parent = baseenv())
  sys.source(user_file, envir = user)
  if (!is.list(user$config)) stop("Configuration must define a list named config.")
  unknown <- setdiff(names(user$config), names(config))
  if (length(unknown)) stop("Unknown config fields: ", paste(unknown, collapse = ", "))
  config <- utils::modifyList(config, user$config, keep.null = TRUE)
}
dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
config$work_dir <- normalizePath(config$work_dir, mustWork = TRUE)
if (is.null(config$input_dir)) config$input_dir <- file.path(config$work_dir, "inputs")
if (length(config$library_paths)) {
  .libPaths(c(config$library_paths, .libPaths()))
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
}
Sys.setenv(IMMUCELLAI2_FIG3_HOME = home,
           IMMUCELLAI2_FIG3_WORK = config$work_dir,
           IMMUCELLAI2_FIG3_INPUT = normalizePath(config$input_dir, mustWork = FALSE),
           IMMUCELLAI2_FIG2_RESULTS = normalizePath(config$fig2_results_dir, mustWork = FALSE),
           IMMUCELLAI2_FIG3_CORES = config$n_cores,
           IMMUCELLAI2_FIG3_REFERENCE = if (is.null(config$reference_file)) "" else config$reference_file,
           IMMUCELLAI2_FIG3_MARKERS = if (is.null(config$immune_gene_file)) "" else config$immune_gene_file,
           IMMUCELLAI2_FIG3_ARCHIVED_TIMING = if (isTRUE(config$use_archived_timing)) "1" else "0")
source(file.path(home, "paths.R"))

run_script <- function(relative) {
  env <- new.env(parent = globalenv())
  env$fig3_requested_datasets <- character(0)
  sys.source(fig3_script(relative), envir = env)
  if (exists("run_status", envir = env, inherits = FALSE) &&
      any(env$run_status$Status == "FAILED")) stop("A method failed: ", relative)
  invisible(NULL)
}
run_stage <- function(stage) {
  message("Figure 3 stage: ", stage)
  switch(stage,
    prepare_pbmc = run_script("pbmc/prepare_newman_monaco_inputs.R"),
    compare_pbmc = run_script("pbmc/run_newman_monaco_7tools_compare.R"),
    compare_xu = run_script("xu/compare_gse107011_result_tpm_7tools_best_immucellai2.R"),
    xu_pairwise = run_script("xu/compare_gse107019_result_tpm_immucellai2_vs_bayesprism.R"),
    evaluate_xu = run_script("xu/recalculate_plot_gse107011_adjusted_7tools_flatvb.R"),
    plot_xu = {
      env <- new.env(parent = globalenv())
      env$out_dir <- file.path(fig3_work_dir(), "combined_4datasets_7tools_pie_plots_final_ordered")
      env$input_files <- data.frame(Dataset = "GSE107011", File = file.path(fig3_work_dir(),
        "GSE107011_adjusted_7tools_flatvb", "GSE107011_adjusted_7tools_per_celltype_metrics.txt"))
      sys.source(file.path(repo, "benchmarks", "fig2",
        "plot_4datasets_7tools_correlation_pie_matrix_final_ordered.R"), envir = env)
    },
    plot_pbmc = run_script("pbmc/plot_newman_monaco_7tools_per_celltype_pies.R"),
    runtime = run_script("runtime/run_runtime_benchmark_1sample_vs_50samples.R"),
    plot_runtime = run_script("runtime/plot_runtime_benchmark_1sample_vs_50samples_logvalue.R"),
    runtime_summary = {
      run_script("runtime/summarize_runtime_benchmark_1sample_vs_50samples.R")
      run_script("runtime/plot_runtime_benchmark_1sample_vs_50samples_logvalue.R")
    },
    summarize = {
      run_script("summary/plot_4datasets_7tools_overall_bar_metrics.R")
      run_script("summary/plot_6datasets_7tools_mean_pearson.R")
    },
    compose = run_script("plots/assemble_current_figure3.R")
  )
}
sequence <- if (args[1] == "all") {
  c("prepare_pbmc", "compare_xu", "xu_pairwise", "evaluate_xu", "compare_pbmc",
    "runtime", "runtime_summary", "plot_xu", "plot_pbmc", "summarize", "compose")
} else if (args[1] == "plots") {
  c("plot_xu", "plot_pbmc", "plot_runtime", "summarize", "compose")
} else args[1]
tryCatch({
  for (stage in sequence) run_stage(stage)
}, finally = {
  writeLines(capture.output(sessionInfo()), file.path(config$work_dir, "sessionInfo.txt"))
  saveRDS(config, file.path(config$work_dir, "run_config.rds"))
})
