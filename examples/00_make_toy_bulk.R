suppressPackageStartupMessages(library(ImmuCellAI2.0))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1] else "example_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reference <- load_immucellai2_reference()
markers <- load_immucellai2_markers()
genes <- intersect(markers, rownames(reference))

states <- c("CD4Tn", "CD8Tem", "Bnaive", "cMo")
if (!all(states %in% colnames(reference))) {
  stop("The packaged reference is missing a state required by the toy example.")
}

truth <- rbind(
  Toy_1 = c(0.55, 0.20, 0.15, 0.10),
  Toy_2 = c(0.20, 0.50, 0.10, 0.20),
  Toy_3 = c(0.10, 0.15, 0.55, 0.20),
  Toy_4 = c(0.15, 0.10, 0.20, 0.55)
)
colnames(truth) <- states

bulk <- reference[genes, states, drop = FALSE] %*% t(truth)
colnames(bulk) <- rownames(truth)

write_matrix <- function(x, file, row_label) {
  out <- data.frame(
    gene = rownames(x),
    x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- row_label
  utils::write.table(
    out,
    file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_matrix(bulk, file.path(out_dir, "toy_bulk_TPM.txt"), "gene")
write_matrix(truth, file.path(out_dir, "toy_truth_fraction.txt"), "sample")

cat("Toy bulk:", normalizePath(file.path(out_dir, "toy_bulk_TPM.txt")), "\n")
cat("Toy truth:", normalizePath(file.path(out_dir, "toy_truth_fraction.txt")), "\n")
cat("Genes:", nrow(bulk), " Samples:", ncol(bulk), "\n")
