args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript examples/04_custom_reference.R bulk_TPM.txt custom_reference_TPM.txt [output_dir] [cores]")
}
output_dir <- if (length(args) >= 3L) args[3] else "custom_reference_results"
n_cores <- if (length(args) >= 4L) as.integer(args[4]) else 1L

fit <- ImmuCellAI2.0::run_immucellai2_custom(
  bulk = args[1], reference = args[2], n.cores = n_cores
)
ImmuCellAI2.0::write_deconvolution_outputs(fit, output_dir)
saveRDS(fit$settings, file.path(output_dir, "settings.rds"))
writeLines(fit$genes.used, file.path(output_dir, "genes_used.txt"))
cat("Samples:", nrow(fit$subtype.fraction), "\n")
cat("Reference populations:", ncol(fit$subtype.fraction), "\n")
cat("Shared genes used:", length(fit$genes.used), "\n")
cat("Results:", normalizePath(output_dir, mustWork = FALSE), "\n")
