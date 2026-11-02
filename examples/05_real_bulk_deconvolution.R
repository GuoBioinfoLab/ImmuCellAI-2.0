# Run from the repository root in R or RStudio.
bulk <- ImmuCellAI2.0::read_expression_matrix("example_data/Example_expression.txt")
print(dim(bulk))
# 36192 genes x 18 samples

fit <- ImmuCellAI2.0::run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 2
)

print(dim(fit$subtype.fraction))
# 18 samples x 53 terminal immune subtypes
print(dim(fit$type.fraction))
print(dim(fit$lineage.fraction))
print(head(fit$subtype.fraction[, 1:6, drop = FALSE]))
print(rowSums(fit$subtype.fraction))

output_dir <- "example_output/real_bulk_ImmuCellAI2"
ImmuCellAI2.0::write_deconvolution_outputs(fit, output_dir)
cat("Results:", normalizePath(output_dir, mustWork = TRUE), "\n")
invisible(fit)
