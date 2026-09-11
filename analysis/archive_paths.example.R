# Copy to analysis/archive_paths.R and use absolute paths with forward slashes.
# LOCAL_R_ROOT contains Fig4/, Fig5/, Fig6/, the reference and marker files.
path_roots <- c(
  LOCAL_R_ROOT = "/path/to/local/R",
  LOCAL_CLUSTER_ROOT = "/path/to/cluster_TCGA",
  LOCAL_PROJECT_ROOT = "/path/to/project",
  LOCAL_LEGACY_PROJECT_ROOT = "/path/to/legacy/project",
  LOCAL_LEGACY_AUX_ROOT = "/path/to/legacy/auxiliary"
)
prepared_dir <- "analysis/local_sources"
