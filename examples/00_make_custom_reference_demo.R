args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "example_data/custom_reference"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Synthetic identifiers and signatures demonstrate the input format, not biology.
set.seed(29)
reference <- matrix(rexp(1500 * 3), nrow = 1500)
dimnames(reference) <- list(sprintf("DemoGene_%04d", 1:1500),
                           c("Population_A", "Population_B", "Population_C"))
for (k in 1:3) {
  genes <- seq.int(500 * (k - 1) + 1, 500 * k)
  reference[genes, k] <- reference[genes, k] + 15
}
reference <- sweep(reference, 2, colSums(reference), "/") * 1e6
truth <- rbind(
  Custom_1 = c(0.70, 0.20, 0.10),
  Custom_2 = c(0.20, 0.30, 0.50),
  Custom_3 = c(0.10, 0.75, 0.15),
  Custom_4 = c(0.35, 0.15, 0.50)
)
colnames(truth) <- colnames(reference)
bulk <- reference %*% t(truth)

write_matrix <- function(x, filename, row_label) {
  tab <- data.frame(id = rownames(x), x, check.names = FALSE)
  names(tab)[1] <- row_label
  utils::write.table(tab, file.path(out_dir, filename), sep = "\t",
                     quote = FALSE, row.names = FALSE)
}
write_matrix(bulk, "bulk_TPM.txt", "gene")
write_matrix(reference, "custom_reference_TPM.txt", "gene")
write_matrix(truth, "truth_fraction.txt", "sample")
cat("Synthetic custom-reference demo:", normalizePath(out_dir), "\n")
