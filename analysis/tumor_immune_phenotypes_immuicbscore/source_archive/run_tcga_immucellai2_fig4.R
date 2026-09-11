options(stringsAsFactors = FALSE)

input_file <- "<LOCAL_R_ROOT>/Fig4/TCGA_TPM_symbol.txt"
reference_file <- "<LOCAL_R_ROOT>/reference_53celltypesTPM20260518.txt"
immune_gene_file <- "<LOCAL_R_ROOT>/MarkerUsedDeconvolution.txt"
out_dir <- "<LOCAL_R_ROOT>/Fig4"
filtered_file <- file.path(out_dir, "TCGA_TPM_symbol_immune_genes.txt")
checkpoint_dir <- file.path(out_dir, "ImmuCellAI2_TCGA_checkpoints")
stable_tmp <- file.path(out_dir, "ImmuCellAI2_TCGA_tmp")

n_cores <- min(8L, parallel::detectCores())
chunk_size <- 400L
set.seed(123)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

suppressPackageStartupMessages({
  library(ImmuCellAI2.0)
  library(data.table)
})

log_file <- file.path(out_dir, "ImmuCellAI2_TCGA_run.log")
log_message <- function(...) {
  text <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(text, "\n")
  cat(text, "\n", file = log_file, append = TRUE)
}

read_gene_list <- function(file) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[\t,; ]+"))
  x <- trimws(x)
  unique(x[nzchar(x)])
}

filter_expression_file <- function(input, output, genes) {
  if (file.exists(output) && file.info(output)$size > 1e6) {
    log_message("Reusing filtered immune-gene matrix: ", output)
    return(invisible(output))
  }

  log_message("Streaming input and retaining ", length(genes), " candidate immune genes...")
  genes_env <- new.env(hash = TRUE, parent = emptyenv())
  for (g in genes) assign(g, TRUE, envir = genes_env)

  con_in <- file(input, open = "r", encoding = "UTF-8")
  con_out <- file(output, open = "w", encoding = "UTF-8")
  on.exit({
    try(close(con_in), silent = TRUE)
    try(close(con_out), silent = TRUE)
  }, add = TRUE)

  header <- readLines(con_in, n = 1L, warn = FALSE)
  writeLines(paste0("Gene\t", header), con_out, sep = "\n")

  retained <- 0L
  scanned <- 0L
  repeat {
    lines <- readLines(con_in, n = 100L, warn = FALSE)
    if (!length(lines)) break
    first_tabs <- regexpr("\t", lines, fixed = TRUE)
    gene_names <- substring(lines, 1L, first_tabs - 1L)
    keep <- vapply(gene_names, exists, logical(1), envir = genes_env, inherits = FALSE)
    if (any(keep)) {
      writeLines(lines[keep], con_out, sep = "\n")
      retained <- retained + sum(keep)
    }
    scanned <- scanned + length(lines)
    if (scanned %% 5000L == 0L) log_message("Scanned genes: ", scanned, "; retained rows: ", retained)
  }
  log_message("Filtering complete. Scanned genes: ", scanned, "; retained rows: ", retained)
  invisible(output)
}

aggregate_duplicate_genes <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  log_message("Aggregating duplicated gene symbols...")
  rowsum(mat, group = rownames(mat), reorder = FALSE)
}

log_message("Reading reference...")
reference <- read_expression_matrix(reference_file)
immune_genes <- intersect(read_gene_list(immune_gene_file), rownames(reference))
if (length(immune_genes) < 100L) stop("Too few immune genes overlap the reference.")

filter_expression_file(input_file, filtered_file, immune_genes)

log_message("Reading filtered TCGA matrix...")
bulk_dt <- fread(filtered_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
gene_names <- as.character(bulk_dt[[1L]])
bulk_dt[[1L]] <- NULL
bulk <- as.matrix(bulk_dt)
storage.mode(bulk) <- "double"
rownames(bulk) <- gene_names
rm(bulk_dt, gene_names)
gc()

bulk <- aggregate_duplicate_genes(bulk)
common <- intersect(rownames(bulk), rownames(reference))
bulk <- bulk[common, , drop = FALSE]
reference <- reference[common, , drop = FALSE]
keep <- rowSums(reference) > 0 & rowSums(bulk) > 0
bulk <- bulk[keep, , drop = FALSE]
reference <- reference[keep, , drop = FALSE]

hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference)
reference <- reference[, hierarchy$state_name, drop = FALSE]

samples <- colnames(bulk)
chunks <- split(seq_along(samples), ceiling(seq_along(samples) / chunk_size))
log_message(
  "Prepared matrix: ", nrow(bulk), " genes x ", ncol(bulk), " samples; ",
  length(chunks), " chunks; ", n_cores, " cores."
)

status <- data.frame(
  Chunk = seq_along(chunks),
  Start = vapply(chunks, min, integer(1)),
  End = vapply(chunks, max, integer(1)),
  Samples = vapply(chunks, length, integer(1)),
  Status = "PENDING",
  Seconds = NA_real_,
  stringsAsFactors = FALSE
)

state_parts <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  idx <- chunks[[i]]
  checkpoint <- file.path(checkpoint_dir, sprintf("chunk_%03d_state_fraction.rds", i))
  if (file.exists(checkpoint)) {
    log_message("Reusing checkpoint chunk ", i, "/", length(chunks))
    state_parts[[i]] <- readRDS(checkpoint)
    status$Status[i] <- "REUSED"
    next
  }

  log_message("Running chunk ", i, "/", length(chunks), " (", length(idx), " samples)...")
  status$Status[i] <- "RUNNING"
  write.table(status, file.path(out_dir, "ImmuCellAI2_TCGA_chunk_status.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  started <- proc.time()[["elapsed"]]

  fit <- deconvolve_bulk_matrix(
    bulk.mat = bulk[, idx, drop = FALSE],
    reference = reference,
    hierarchy = hierarchy,
    hierarchy.mode = "tcell",
    inference.method = "vb",
    add.unknown = FALSE,
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    n.cores = n_cores,
    seed = 123 + i
  )

  state_parts[[i]] <- as.matrix(fit$state.fraction)
  saveRDS(state_parts[[i]], checkpoint, compress = FALSE)
  status$Status[i] <- "OK"
  status$Seconds[i] <- proc.time()[["elapsed"]] - started
  write.table(status, file.path(out_dir, "ImmuCellAI2_TCGA_chunk_status.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_message("Finished chunk ", i, " in ", round(status$Seconds[i], 1), " seconds.")
  rm(fit)
  gc()
}

state_fraction_observed <- do.call(rbind, state_parts)
state_names <- colnames(state_fraction_observed)
state_fraction_samples_by_states <- matrix(
  NA_real_,
  nrow = length(samples),
  ncol = length(state_names),
  dimnames = list(samples, state_names)
)
matched_samples <- intersect(samples, rownames(state_fraction_observed))
state_fraction_samples_by_states[matched_samples, ] <-
  state_fraction_observed[matched_samples, state_names, drop = FALSE]

missing_samples <- setdiff(samples, rownames(state_fraction_observed))
write.table(
  data.frame(
    Sample = missing_samples,
    Reason = "zero expression across all retained immune genes",
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "TCGA_ImmuCellAI2_unprocessed_samples.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

major_groups <- split(hierarchy$state_name, hierarchy$major_lineage)
major_fraction_samples_by_lineage <- sapply(major_groups, function(states) {
  rowSums(state_fraction_samples_by_states[, states, drop = FALSE], na.rm = FALSE)
})
major_fraction_samples_by_lineage <- as.matrix(major_fraction_samples_by_lineage)

state_fraction <- t(state_fraction_samples_by_states)
major_fraction <- t(major_fraction_samples_by_lineage)

write.table(
  state_fraction_samples_by_states,
  file.path(out_dir, "TCGA_ImmuCellAI2_state_fraction.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(
  major_fraction_samples_by_lineage,
  file.path(out_dir, "TCGA_ImmuCellAI2_major_lineage_fraction.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(
  hierarchy,
  file.path(out_dir, "TCGA_ImmuCellAI2_hierarchy.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

run_info <- data.frame(
  Package = "ImmuCellAI2.0",
  Version = as.character(packageVersion("ImmuCellAI2.0")),
  Samples = ncol(bulk),
  GenesUsed = nrow(bulk),
  States = nrow(state_fraction),
  SamplesProcessed = length(matched_samples),
  SamplesUnprocessed = length(missing_samples),
  HierarchyMode = "tcell",
  InferenceMethod = "vb",
  AddUnknown = FALSE,
  Cores = n_cores,
  ChunkSize = chunk_size,
  NIter = 50,
  VBTolerance = 1e-6,
  stringsAsFactors = FALSE
)
write.table(run_info, file.path(out_dir, "TCGA_ImmuCellAI2_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(
    state.fraction = state_fraction,
    major.fraction = major_fraction,
    hierarchy = hierarchy,
    run.info = run_info
  ),
  file.path(out_dir, "TCGA_ImmuCellAI2_results.rds"),
  compress = FALSE
)

log_message("All outputs completed.")
print(run_info)
