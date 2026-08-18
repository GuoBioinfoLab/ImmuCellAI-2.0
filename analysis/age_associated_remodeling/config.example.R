repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_root, "analysis", "age_associated_remodeling", "data")
output_dir <- file.path(repo_root, "analysis", "age_associated_remodeling", "results")

expression_files <- c(
  file.path(data_dir, "age_bulk.txt"),
  file.path(data_dir, "age_bulk2.txt")
)
metadata_files <- c(
  file.path(data_dir, "age_sample_info.csv"),
  file.path(data_dir, "healthy_info2_clean.csv")
)
sample_column_candidates <- c("sample", "Sample", "sample_id", "SampleID", "Run", "ID")
age_column_candidates <- c("age", "Age", "age_years", "Age_years")
n_cores <- 8L
