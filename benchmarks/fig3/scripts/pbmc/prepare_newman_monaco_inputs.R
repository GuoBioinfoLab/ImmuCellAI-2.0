# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

root <- file.path(fig3_work_dir(), "Fig2/Newman_Monaco_7tools_validation")

write_standardized <- function(dataset, bulk_file, truth_file, sample_transform = identity) {
  out_dir <- file.path(root, dataset, "standardized")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bulk <- read.delim(bulk_file, check.names = FALSE)
  rownames(bulk) <- bulk[[1]]
  bulk[[1]] <- NULL
  colnames(bulk) <- sample_transform(colnames(bulk))

  truth <- read.csv(truth_file, check.names = FALSE)
  rownames(truth) <- as.character(truth[[1]])
  truth[[1]] <- NULL

  samples <- intersect(colnames(bulk), rownames(truth))
  if (length(samples) < 3L) stop("Sample alignment failed for ", dataset)
  bulk <- bulk[, samples, drop = FALSE]
  truth <- truth[samples, , drop = FALSE]

  write.csv(data.frame(Gene = rownames(bulk), bulk, check.names = FALSE),
            file.path(out_dir, "expression.csv"), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(Sample = rownames(truth), truth, check.names = FALSE),
            file.path(out_dir, "true_proportions.csv"), row.names = FALSE, quote = FALSE)
}

write_standardized(
  "RNA_Sieve_newman_pbmcs",
  file.path(root, "RNA_Sieve_newman_pbmcs", "newman_bulk_data.txt"),
  file.path(root, "RNA_Sieve_newman_pbmcs", "newman_true_props.csv")
)

write_standardized(
  "RNA_Sieve_monaco_pbmcs",
  file.path(root, "RNA_Sieve_monaco_pbmcs", "monaco_bulk_data.txt"),
  file.path(root, "RNA_Sieve_monaco_pbmcs", "monaco_true_props.csv"),
  function(x) sub("_PBMC$", "", x)
)
