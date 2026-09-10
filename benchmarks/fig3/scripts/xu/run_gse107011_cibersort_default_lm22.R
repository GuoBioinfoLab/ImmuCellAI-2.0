# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(CIBERSORT))

bulk_file <- file.path(fig3_input_dir(), "result_tpm.csv")
sample_id_file <- file.path(fig3_input_dir(), "GSE107019_sampleID.csv")
out_dir <- file.path(fig3_work_dir(), "GSE107011_result_tpm_7tools_best_immucellai2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_bulk <- function(file) {
  x <- read.delim(file, row.names = 1, check.names = FALSE)
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

bulk <- read_bulk(bulk_file)
sample_map <- read.csv(sample_id_file, header = FALSE, check.names = FALSE)
bulk_samples <- trimws(as.character(sample_map[[1]]))
bulk_samples <- bulk_samples[bulk_samples %in% colnames(bulk)]
bulk <- bulk[, bulk_samples, drop = FALSE]

data(LM22, package = "CIBERSORT")
genes <- intersect(rownames(LM22), rownames(bulk))
message("Running CIBERSORT default LM22: ", length(genes), "/547 genes; ",
        ncol(bulk), " samples")

# Package defaults: maxSize=500, perm=0, QN=TRUE.
fit <- CIBERSORT::cibersort(
  sig_matrix = LM22[genes, , drop = FALSE],
  mixture_file = bulk[genes, , drop = FALSE]
)
fit <- as.matrix(fit)
fraction_cols <- intersect(colnames(LM22), colnames(fit))
fraction <- fit[bulk_samples, fraction_cols, drop = FALSE]

write.table(
  fraction,
  file.path(out_dir, "CIBERSORT_default_LM22_fraction.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)
writeLines(
  c(
    "CIBERSORT default LM22",
    paste0("LM22 genes available: ", length(genes), "/547"),
    paste0("Aligned bulk samples: ", ncol(bulk)),
    "Parameters: package defaults maxSize=500, perm=0, QN=TRUE"
  ),
  file.path(out_dir, "CIBERSORT_default_LM22_run_info.txt")
)

cat("Saved:", file.path(out_dir, "CIBERSORT_default_LM22_fraction.txt"), "\n")
