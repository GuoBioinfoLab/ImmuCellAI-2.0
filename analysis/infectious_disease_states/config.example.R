repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_root, "analysis", "infectious_disease_states", "data")
output_dir <- file.path(repo_root, "analysis", "infectious_disease_states", "results")

expression_file <- file.path(data_dir, "exp.txt")
tb_metadata_file <- file.path(data_dir, "PRJNA_combined.csv")
hiv_metadata_file <- file.path(data_dir, "PRJNA683803_combined.csv")
hiv_longitudinal_file <- file.path(data_dir, "PRJNA683803.xlsx")
fraction_file <- file.path(output_dir, "deconvolution", "state_fraction.txt")
n_cores <- 8L
