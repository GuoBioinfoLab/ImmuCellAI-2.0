options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
workspace_dir <- "<LOCAL_R_ROOT>"

clinical_file <- file.path(fig4_dir, "clinical.tsv")
tcga_tpm_file <- file.path(fig4_dir, "TCGA_TPM_symbol.txt")
reference_file <- file.path(workspace_dir, "reference_53celltypesTPM20260518.txt")
marker_file <- file.path(workspace_dir, "MarkerUsedDeconvolution.txt")

out_dir <- file.path(fig4_dir, "TCGA_BRCA_ImmuCellAI2_marker9000")
checkpoint_dir <- file.path(out_dir, "checkpoints")
tmp_dir <- file.path(out_dir, "tmp")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = tmp_dir, TEMP = tmp_dir, TMPDIR = tmp_dir)

n_cores <- min(8L, parallel::detectCores(logical = TRUE))
chunk_size <- 250L

suppressPackageStartupMessages({
  library(data.table)
  library(ImmuCellAI2.0)
})

log_file <- file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_run.log")
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

read_gene_tokens <- function(file) {
  txt <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(txt, "[\t,;\"' ]+"))
  x <- trimws(x)
  x <- unique(x[nzchar(x)])
  x <- x[!grepl("^[0-9]+$", x)]
  x <- x[!x %in% c("Gene", "Genes", "gene", "genes", "X", "x", "NA")]
  x
}

read_matrix_gene_first <- function(file) {
  con <- file(file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  second <- readLines(con, n = 1L, warn = FALSE)
  first_fields <- strsplit(first, "\t", fixed = TRUE)[[1L]]
  second_fields <- strsplit(second, "\t", fixed = TRUE)[[1L]]
  if (length(first_fields) == length(second_fields) - 1L) {
    dt <- fread(
      file,
      sep = "\t",
      header = FALSE,
      skip = 1L,
      col.names = c("Gene", first_fields),
      check.names = FALSE,
      data.table = FALSE
    )
  } else {
    dt <- fread(file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  }
  gene <- as.character(dt[[1L]])
  dt[[1L]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "double"
  rownames(mat) <- gene
  mat[!is.finite(mat)] <- 0
  mat
}

aggregate_duplicate_genes <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  rowsum(mat, group = rownames(mat), reorder = FALSE)
}

sample_type_label <- function(code) {
  lookup <- c(
    "01" = "Primary Tumor",
    "02" = "Recurrent Solid Tumor",
    "03" = "Primary Blood Derived Cancer",
    "05" = "Additional New Primary",
    "06" = "Metastatic",
    "10" = "Blood Derived Normal",
    "11" = "Solid Tissue Normal"
  )
  out <- unname(lookup[code])
  out[is.na(out)] <- "Other"
  out
}

log_message("Reading clinical file...")
clinical <- fread(clinical_file, sep = "\t", data.table = FALSE)
brca_cases <- unique(clinical$case_submitter_id[clinical$project_id == "TCGA-BRCA"])
brca_cases <- brca_cases[nzchar(brca_cases)]
if (length(brca_cases) == 0L) stop("No TCGA-BRCA cases found in clinical.tsv")

log_message("Reading TCGA TPM header...")
all_samples <- scan(tcga_tpm_file, what = "", nlines = 1L, sep = "\t", quiet = TRUE)
brca_samples <- all_samples[substr(all_samples, 1L, 12L) %in% brca_cases]
if (length(brca_samples) == 0L) stop("No BRCA samples matched between clinical.tsv and TCGA_TPM_symbol.txt")

sample_info <- data.frame(
  Sample = brca_samples,
  CaseSubmitterID = substr(brca_samples, 1L, 12L),
  ProjectID = "TCGA-BRCA",
  SampleTypeCode = substr(brca_samples, 14L, 15L),
  SampleType = sample_type_label(substr(brca_samples, 14L, 15L)),
  stringsAsFactors = FALSE
)
sample_info <- unique(sample_info)
write.table(sample_info, file.path(out_dir, "TCGA_BRCA_sample_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

log_message("BRCA cases in clinical: ", length(brca_cases))
log_message("Matched BRCA expression samples: ", length(brca_samples))
log_message("Sample type counts: ", paste(names(table(sample_info$SampleTypeCode)), table(sample_info$SampleTypeCode), sep = "=", collapse = "; "))

select_cols <- c(1L, match(brca_samples, all_samples) + 1L)
brca_expr_file <- file.path(out_dir, "TCGA_BRCA_TPM_symbol.txt")
if (!file.exists(brca_expr_file)) {
  log_message("Extracting BRCA expression matrix from full TCGA TPM file...")
  brca_dt <- fread(
    tcga_tpm_file,
    sep = "\t",
    header = FALSE,
    skip = 1L,
    select = select_cols,
    col.names = c("Gene", brca_samples),
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE
  )
  fwrite(brca_dt, brca_expr_file, sep = "\t", quote = FALSE)
  rm(brca_dt)
  gc()
} else {
  log_message("Reusing BRCA expression matrix: ", brca_expr_file)
}

log_message("Reading reference and BRCA matrix...")
reference <- read_matrix_gene_first(reference_file)
bulk <- read_matrix_gene_first(brca_expr_file)
reference <- aggregate_duplicate_genes(reference)
bulk <- aggregate_duplicate_genes(bulk)

marker_genes <- read_gene_tokens(marker_file)
common <- Reduce(intersect, list(rownames(reference), rownames(bulk), marker_genes))
if (length(common) < 100L) stop("Too few marker genes overlap BRCA bulk and reference: ", length(common))

bulk <- bulk[common, , drop = FALSE]
reference <- reference[common, , drop = FALSE]
keep <- rowSums(reference) > 0 & rowSums(bulk) > 0
bulk <- bulk[keep, , drop = FALSE]
reference <- reference[keep, , drop = FALSE]

hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference)
reference <- reference[, hierarchy$state_name, drop = FALSE]

marker_filtered_file <- file.path(out_dir, "TCGA_BRCA_TPM_symbol_marker5000_filtered.txt")
if (!file.exists(marker_filtered_file)) {
  write.table(
    data.frame(Gene = rownames(bulk), bulk, check.names = FALSE),
    marker_filtered_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

samples <- colnames(bulk)
chunks <- split(seq_along(samples), ceiling(seq_along(samples) / chunk_size))
log_message("Prepared matrix: ", nrow(bulk), " genes x ", ncol(bulk), " samples; chunks=", length(chunks), "; cores=", n_cores)

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
  write.table(status, file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_chunk_status.txt"),
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
    seed = 123 + i,
    verbose = FALSE
  )

  part <- as.matrix(fit$state.fraction)
  colnames(part) <- hierarchy$state_name
  state_parts[[i]] <- part
  saveRDS(part, checkpoint, compress = FALSE)

  status$Status[i] <- "OK"
  status$Seconds[i] <- proc.time()[["elapsed"]] - started
  write.table(status, file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_chunk_status.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_message("Finished chunk ", i, " in ", round(status$Seconds[i], 1), " seconds.")
  rm(fit, part)
  gc()
}

state_fraction_samples_by_states <- do.call(rbind, state_parts)
state_fraction_samples_by_states <- state_fraction_samples_by_states[samples, hierarchy$state_name, drop = FALSE]

major_groups <- split(hierarchy$state_name, hierarchy$major_lineage)
major_fraction_samples_by_lineage <- sapply(major_groups, function(states) {
  rowSums(state_fraction_samples_by_states[, states, drop = FALSE], na.rm = FALSE)
})
major_fraction_samples_by_lineage <- as.matrix(major_fraction_samples_by_lineage)

write.table(
  state_fraction_samples_by_states,
  file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  t(state_fraction_samples_by_states),
  file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_state_fraction_celltype_by_sample.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  major_fraction_samples_by_lineage,
  file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_major_fraction_sample_by_lineage.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(hierarchy, file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_hierarchy.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

run_info <- data.frame(
  Package = "ImmuCellAI2.0",
  Version = as.character(packageVersion("ImmuCellAI2.0")),
  ProjectID = "TCGA-BRCA",
  ClinicalCases = length(brca_cases),
  Samples = ncol(bulk),
  GenesUsed = nrow(bulk),
  States = ncol(state_fraction_samples_by_states),
  HierarchyMode = "tcell",
  InferenceMethod = "vb",
  AddUnknown = FALSE,
  GeneSet = "Immune 5000+ MarkerUsedDeconvolution",
  Cores = n_cores,
  ChunkSize = chunk_size,
  NIter = 50,
  VBTolerance = 1e-6,
  stringsAsFactors = FALSE
)
write.table(run_info, file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(
    state.fraction = state_fraction_samples_by_states,
    major.fraction = major_fraction_samples_by_lineage,
    sample.info = sample_info,
    hierarchy = hierarchy,
    run.info = run_info
  ),
  file.path(out_dir, "TCGA_BRCA_ImmuCellAI2_marker5000_results.rds"),
  compress = FALSE
)

log_message("All outputs completed.")
print(run_info)
