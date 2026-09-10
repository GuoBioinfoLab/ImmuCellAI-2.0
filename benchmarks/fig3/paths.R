fig3_home <- function() normalizePath(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), mustWork = TRUE)
fig3_work_dir <- function() Sys.getenv("IMMUCELLAI2_FIG3_WORK", "benchmarks/output/fig3")
fig3_input_dir <- function() Sys.getenv("IMMUCELLAI2_FIG3_INPUT", file.path(fig3_work_dir(), "inputs"))
fig3_script <- function(path) file.path(fig3_home(), "scripts", path)
fig3_tumor_metrics <- function(dataset) {
  file.path(Sys.getenv("IMMUCELLAI2_FIG2_RESULTS", "benchmarks/output/fig2"),
            dataset, "all_methods_per_celltype_metrics.txt")
}
fig3_reference <- function() {
  path <- Sys.getenv("IMMUCELLAI2_FIG3_REFERENCE", "")
  if (!nzchar(path)) path <- system.file("extdata", "reference_53celltypesTPM20260518.txt", package = "ImmuCellAI2.0")
  if (!file.exists(path)) stop("Reference atlas missing; install ImmuCellAI2.0 or configure reference_file.")
  path
}
fig3_markers <- function() {
  path <- Sys.getenv("IMMUCELLAI2_FIG3_MARKERS", "")
  if (!nzchar(path)) path <- system.file("extdata", "MarkerUsedDeconvolution_5510.txt", package = "ImmuCellAI2.0")
  if (!file.exists(path)) stop("Marker panel missing; install ImmuCellAI2.0 or configure immune_gene_file.")
  path
}
fig3_cores <- function() {
  value <- as.integer(Sys.getenv("IMMUCELLAI2_FIG3_CORES", "8"))
  if (length(value) != 1L || !is.finite(value) || value < 1L) stop("Invalid core count.")
  value
}
fig3_run_prism <- function(...) {
  # Some BayesPrism versions remove the session temporary directory on exit.
  on.exit(invisible(tempdir(check = TRUE)), add = TRUE)
  BayesPrism::run.prism(...)
}
