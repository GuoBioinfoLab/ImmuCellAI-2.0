args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("ImmuCellAI2.0", "data.table"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_files <- vapply(expression_files, assert_file, character(1), label = "Age expression matrix")
matrices <- lapply(expression_files, ImmuCellAI2.0::read_expression_matrix)
common_genes <- Reduce(intersect, lapply(matrices, rownames))
if (length(common_genes) < 1000L) stop("Too few shared genes across age expression matrices: ", length(common_genes))
bulk <- do.call(cbind, lapply(matrices, function(x) x[common_genes, , drop = FALSE]))
if (anyDuplicated(colnames(bulk))) stop("Duplicate sample identifiers across age expression matrices.")

metadata <- do.call(rbind, lapply(metadata_files[file.exists(metadata_files)], function(path) {
  x <- read_table_auto(path)
  sample_col <- intersect(sample_column_candidates, names(x))[1L]
  age_col <- intersect(age_column_candidates, names(x))[1L]
  if (is.na(sample_col) || is.na(age_col)) {
    stop("Metadata must contain a recognized sample and age column: ", path)
  }
  data.frame(sample = as.character(x[[sample_col]]), age = as.numeric(x[[age_col]]))
}))
metadata <- metadata[!duplicated(metadata$sample) & is.finite(metadata$age), ]
keep <- intersect(colnames(bulk), metadata$sample)
if (!length(keep)) stop("No sample names overlap expression and age metadata.")
bulk <- bulk[, keep, drop = FALSE]
metadata <- metadata[match(keep, metadata$sample), ]

combined_file <- file.path(output_dir, "healthy_age_combined_TPM.tsv")
data.table::fwrite(
  data.frame(gene = rownames(bulk), bulk, check.names = FALSE),
  combined_file, sep = "\t", quote = FALSE
)
data.table::fwrite(metadata, file.path(output_dir, "healthy_age_aligned_metadata.tsv"), sep = "\t")
run_standard_deconvolution(combined_file, file.path(output_dir, "deconvolution"), n_cores)
