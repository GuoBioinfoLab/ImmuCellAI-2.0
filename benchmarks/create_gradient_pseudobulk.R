suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage: Rscript benchmarks/create_gradient_pseudobulk.R",
      "celltype_mean_TPM.txt output_dir [target_cell_file] [seed]"
    )
  )
}

reference_file <- args[1]
output_dir <- args[2]
target_file <- if (length(args) >= 3L && nzchar(args[3])) args[3] else NA_character_
seed <- if (length(args) >= 4L) as.integer(args[4]) else 123L
if (!is.finite(seed)) seed <- 123L
set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_gene_matrix <- function(file) {
  x <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    row.names = NULL
  )
  genes <- as.character(x[[1]])
  value <- x[-1]
  value[] <- lapply(value, function(column) suppressWarnings(as.numeric(column)))
  mat <- as.matrix(value)
  rownames(mat) <- genes
  colnames(mat) <- names(value)
  if (any(!is.finite(mat)) || any(mat < 0)) {
    stop("Reference values must be finite and non-negative.")
  }
  if (anyDuplicated(genes)) {
    stop("Reference gene symbols must be unique.")
  }
  mat
}

write_gene_matrix <- function(mat, file) {
  out <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
  fwrite(out, file, sep = "\t", quote = FALSE)
}

reference <- read_gene_matrix(reference_file)
if (ncol(reference) < 2L) stop("At least two cell types are required.")

targets <- colnames(reference)
if (!is.na(target_file)) {
  targets <- trimws(readLines(target_file, warn = FALSE))
  targets <- targets[nzchar(targets)]
}
missing_targets <- setdiff(targets, colnames(reference))
if (length(missing_targets)) {
  stop("Targets absent from reference: ", paste(missing_targets, collapse = ", "))
}

gradient <- seq(0.005, 0.2, by = 0.005)
n_samples <- length(targets) * length(gradient)
weights <- matrix(
  0,
  nrow = n_samples,
  ncol = ncol(reference),
  dimnames = list(NULL, colnames(reference))
)
design <- vector("list", n_samples)
sample_ids <- character(n_samples)
row_id <- 1L

for (target in targets) {
  for (step_id in seq_along(gradient)) {
    target_fraction <- gradient[step_id]
    background <- rexp(ncol(reference), rate = 1)
    names(background) <- colnames(reference)
    background[target] <- 0
    background <- background / sum(background) * (1 - target_fraction)

    sample_id <- paste0(target, "_", step_id)
    weights[row_id, ] <- background
    weights[row_id, target] <- target_fraction
    sample_ids[row_id] <- sample_id
    design[[row_id]] <- data.frame(
      Sample = sample_id,
      TargetCell = target,
      Step = step_id,
      TargetProportion = target_fraction,
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
  }
}

rownames(weights) <- sample_ids
bulk <- reference %*% t(weights)
truth_long <- rbindlist(design)
truth_wide <- data.frame(sample = rownames(weights), weights, check.names = FALSE)

write_gene_matrix(bulk, file.path(output_dir, "pseudo_bulk_TPM.txt"))
fwrite(truth_wide, file.path(output_dir, "truth_fraction.txt"), sep = "\t")
fwrite(truth_long, file.path(output_dir, "truth_long.txt"), sep = "\t")
fwrite(
  data.table(
    Seed = seed,
    N_genes = nrow(reference),
    N_cell_types = ncol(reference),
    N_targets = length(targets),
    N_samples = n_samples,
    Gradient = "seq(0.005, 0.2, by = 0.005)"
  ),
  file.path(output_dir, "simulation_settings.txt"),
  sep = "\t"
)

cat("Pseudo-bulk samples:", n_samples, "\n")
cat("Targets:", length(targets), "\n")
cat("Seed:", seed, "\n")
