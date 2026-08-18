suppressPackageStartupMessages(library(ImmuCellAI2.0))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript examples/01_standard_deconvolution.R bulk_TPM.txt [output_dir] [n_cores]")
}

bulk_file <- args[1]
output_dir <- if (length(args) >= 2L) args[2] else "ImmuCellAI2_results"
n_cores <- if (length(args) >= 3L) as.integer(args[3]) else 1L

bulk <- read_expression_matrix(bulk_file)
fit <- run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = n_cores,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-6,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  seed = 123
)

write_deconvolution_outputs(fit, output_dir)
saveRDS(fit$settings, file.path(output_dir, "settings.rds"))
writeLines(fit$genes.used, file.path(output_dir, "genes_used.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

cat("Samples:", nrow(fit$state.fraction), "\n")
cat("Cell states:", ncol(fit$state.fraction), "\n")
cat("Marker genes used:", length(fit$genes.used), "\n")
cat("Results:", normalizePath(output_dir, mustWork = FALSE), "\n")
