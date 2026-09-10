# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

# Compare ImmuCellAI 2.0 and BayesPrism against flow-cytometry
# proportions for GSE107019 / ImmuCellAI2-New.
#
# Inputs:
#   bulk TPM matrix:
#   true flow-cytometry ratios:
#
# Outputs are written to:
#
# Notes:
#   1. result_tpm.csv has a non-standard first line: the sample header has no
#      leading "Gene" field. The reader below handles this format.
#   2. GSE107019_celltypeRatio.csv has a description line followed by the
#      real header line, so it is read with skip = 1.
#   3. If bulk sample IDs and flow-cytometry sample names do not overlap but
#      the sample counts are the same, the script aligns samples by order and
#      writes the order mapping to sample_alignment.txt.

options(stringsAsFactors = FALSE)

bulk_file <- file.path(fig3_input_dir(), "result_tpm.csv")
truth_file <- file.path(fig3_input_dir(), "GSE107019_celltypeRatio.csv")
sample_id_file <- file.path(fig3_input_dir(), "GSE107019_sampleID.csv")
reference_file <- fig3_reference()
out_dir <- file.path(fig3_input_dir(), "GSE107019_compare_ImmuCellAI2_BayesPrism")

n_cores <- fig3_cores()
chain_length <- 600
burn_in <- 500

local_lib <- file.path(fig3_work_dir(), "Rlib_test")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ImmuCellAI2.0", quietly = TRUE)) {
  stop("Install the current ImmuCellAI2.0 package before running this script.")
}
if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  stop(
    "BayesPrism is not installed in this R library. Please install/load BayesPrism, ",
    "then rerun this script. Current .libPaths(): ", paste(.libPaths(), collapse = "; ")
  )
}

library(ImmuCellAI2.0)
library(BayesPrism)

detect_sep <- function(file) {
  x <- readLines(file, n = 1, warn = FALSE)
  if (grepl("\t", x)) "\t" else ","
}

read_bulk_result_tpm <- function(file) {
  sep <- detect_sep(file)
  header <- strsplit(readLines(file, n = 1, warn = FALSE), sep, fixed = TRUE)[[1]]
  header <- trimws(header)
  header <- header[nzchar(header)]

  dat <- read.table(
    file,
    sep = sep,
    header = FALSE,
    skip = 1,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    fill = TRUE
  )
  genes <- as.character(dat[[1]])
  mat <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"

  if (length(header) == ncol(mat)) {
    colnames(mat) <- header
  } else {
    colnames(mat) <- paste0("Sample", seq_len(ncol(mat)))
    warning("Bulk header length did not match matrix columns; using Sample1..SampleN.")
  }
  rownames(mat) <- make.unique(genes)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

read_reference_matrix <- function(file) {
  sep <- detect_sep(file)
  dat <- read.table(
    file,
    sep = sep,
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    row.names = 1
  )
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

read_truth_matrix <- function(file) {
  dat <- read.csv(
    file,
    skip = 1,
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
  sample_col <- grep("^Sample", colnames(dat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sample_col)) sample_col <- colnames(dat)[1]
  sample_names <- as.character(dat[[sample_col]])
  dat <- dat[, setdiff(colnames(dat), sample_col), drop = FALSE]
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- sample_names
  mat[is.na(mat)] <- NA_real_

  # Flow cytometry file stores percentages. Convert to fractions when needed.
  if (max(mat, na.rm = TRUE) > 1.5) mat <- mat / 100
  mat
}

read_sample_id_map <- function(file) {
  if (!file.exists(file)) return(NULL)
  dat <- read.csv(
    file,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  if (ncol(dat) < 2) {
    stop("sample_id_file must contain at least two columns: bulk sample ID and truth sample name.")
  }
  out <- data.frame(
    BulkSample = trimws(as.character(dat[[1]])),
    TruthSample = trimws(as.character(dat[[2]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$BulkSample) & nzchar(out$TruthSample), , drop = FALSE]
  out <- out[!duplicated(out$BulkSample), , drop = FALSE]
  out
}

collapse_duplicate_genes <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  rows <- split(seq_len(nrow(mat)), rownames(mat))
  out <- do.call(rbind, lapply(rows, function(ii) {
    if (length(ii) == 1) mat[ii, ] else colMeans(mat[ii, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(out) <- names(rows)
  out
}

ensure_samples_by_states <- function(pred, sample_names, state_names) {
  pred <- as.matrix(pred)
  if (all(sample_names %in% rownames(pred))) {
    pred <- pred[sample_names, , drop = FALSE]
  } else if (all(sample_names %in% colnames(pred))) {
    pred <- t(pred[, sample_names, drop = FALSE])
  } else if (all(state_names %in% rownames(pred))) {
    pred <- t(pred)
  }
  pred
}

aggregate_prediction <- function(pred, states) {
  states <- intersect(states, colnames(pred))
  if (length(states) == 0) return(rep(NA_real_, nrow(pred)))
  rowSums(pred[, states, drop = FALSE], na.rm = TRUE)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

metric_one <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) == 0) {
    return(c(Pearson = NA, Spearman = NA, RMSE = NA, MAE = NA, Slope = NA, Intercept = NA, N = 0))
  }
  fit <- if (sum(ok) >= 3 && sd(truth[ok]) > 0 && sd(pred[ok]) > 0) {
    coef(lm(pred[ok] ~ truth[ok]))
  } else {
    c(`(Intercept)` = NA_real_, truth = NA_real_)
  }
  c(
    Pearson = safe_cor(pred, truth, "pearson"),
    Spearman = safe_cor(pred, truth, "spearman"),
    RMSE = sqrt(mean((pred[ok] - truth[ok])^2)),
    MAE = mean(abs(pred[ok] - truth[ok])),
    Slope = unname(fit[2]),
    Intercept = unname(fit[1]),
    N = sum(ok)
  )
}

evaluate_against_truth <- function(pred, truth, target_map, method_name) {
  out <- lapply(names(target_map), function(cell) {
    if (!cell %in% colnames(truth)) return(NULL)
    pred_vec <- aggregate_prediction(pred, target_map[[cell]])
    truth_vec <- truth[, cell]
    m <- metric_one(pred_vec, truth_vec)
    data.frame(
      Method = method_name,
      CellType = cell,
      PredictedStates = paste(intersect(target_map[[cell]], colnames(pred)), collapse = "+"),
      Pearson = m["Pearson"],
      Spearman = m["Spearman"],
      RMSE = m["RMSE"],
      MAE = m["MAE"],
      Slope = m["Slope"],
      Intercept = m["Intercept"],
      N_samples = m["N"],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

# Mapping from flow-cytometry columns to ImmuCellAI2/BayesPrism reference states.
# You can edit this list if you want a stricter or broader biological mapping.
target_map <- list(
  "B Naive" = "Bnaive",
  "B NSM" = "FOB",
  "B SM" = "MBC",
  "Plasmablasts" = "PB",
  "T CD4 Naive" = "CD4Tn",
  "Tregs" = "Treg",
  "Tfh" = "Tfh",
  "Th1" = "Th1",
  "Th1/Th17" = "Th1/17",
  "Th17" = "Th17",
  "Th2" = "Th2",
  "T CD8 Naive" = "CD8Tn",
  "T CD8 TE" = "CD8Temra",
  "MAIT" = "MAIT",
  "pDCs" = "pDC",
  "mDCs" = c("cDC1", "cDC2"),
  "Monocytes C" = "cMo",
  "Monocytes I" = "intMo",
  "Neutrophils LD" = "Neutrophil",
  "Basophils LD" = "Basophil",
  "Tfh+Th1-17" = c("Tfh", "Th1", "Th1/17", "Th17"),
  "T gd" = "gdT",
  "Monocytes NC+I" = c("ncMo", "intMo"),
  "Granulocytes LD" = c("Neutrophil", "Basophil", "Eosinophil"),
  "Tfh+Th" = c("Tfh", "Th1", "Th1/17", "Th17", "Th2"),
  "T CD8 Memory" = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
  "DCs" = c("cDC1", "cDC2", "pDC", "monoDC"),
  "Monocytes C+I" = c("cMo", "intMo"),
  "T CD4 Memory" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
  "Monocytes" = c("cMo", "intMo", "ncMo"),
  "T CD4" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  "Myeloid Phagocytes" = c("cMo", "intMo", "ncMo", "M0", "M1", "M2", "TAM", "cDC1", "cDC2", "pDC", "monoDC"),
  "T Naive" = c("CD4Tn", "CD8Tn"),
  "T Memory" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
  "Lymphoid" = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC", "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg", "MAIT", "gdT", "NKT", "cNK", "NKreg", "ILC1", "ILC2", "ILC3"),
  "Myeloid" = c("cMo", "intMo", "ncMo", "M0", "M1", "M2", "TAM", "cDC1", "cDC2", "pDC", "monoDC", "Neutrophil", "Basophil", "Eosinophil", "Mast cell", "MDSC"),
  "Adaptive" = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC", "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  "Innate" = c("MAIT", "gdT", "NKT", "cNK", "NKreg", "ILC1", "ILC2", "ILC3", "cMo", "intMo", "ncMo", "M0", "M1", "M2", "TAM", "cDC1", "cDC2", "pDC", "monoDC", "Neutrophil", "Basophil", "Eosinophil", "Mast cell", "MDSC")
)

message("Reading input matrices...")
bulk <- read_bulk_result_tpm(bulk_file)
reference <- read_reference_matrix(reference_file)
truth <- read_truth_matrix(truth_file)
bulk <- collapse_duplicate_genes(bulk)
reference <- collapse_duplicate_genes(reference)

message("Bulk samples: ", ncol(bulk), "; truth samples: ", nrow(truth))

sample_map <- read_sample_id_map(sample_id_file)
if (!is.null(sample_map)) {
  sample_map <- sample_map[sample_map$BulkSample %in% colnames(bulk) &
                             sample_map$TruthSample %in% rownames(truth), , drop = FALSE]
  if (nrow(sample_map) < 3) {
    stop(
      "Sample ID map was found, but fewer than 3 mapped samples overlap bulk and truth. ",
      "Please check: ", sample_id_file
    )
  }
  sample_map <- sample_map[match(intersect(colnames(bulk), sample_map$BulkSample), sample_map$BulkSample), , drop = FALSE]
  bulk <- bulk[, sample_map$BulkSample, drop = FALSE]
  truth <- truth[sample_map$TruthSample, , drop = FALSE]
  rownames(truth) <- sample_map$BulkSample
  alignment <- sample_map
  message("Samples aligned by sample_id_file: ", nrow(alignment))
} else {
  common_samples <- intersect(colnames(bulk), rownames(truth))
  if (length(common_samples) >= 3) {
  bulk <- bulk[, common_samples, drop = FALSE]
  truth <- truth[common_samples, , drop = FALSE]
  alignment <- data.frame(BulkSample = common_samples, TruthSample = common_samples)
  } else if (ncol(bulk) == nrow(truth)) {
  alignment <- data.frame(BulkSample = colnames(bulk), TruthSample = rownames(truth))
  rownames(truth) <- colnames(bulk)
  warning("No direct sample-name overlap. Samples were aligned by column/row order.")
  } else {
  stop("Cannot align samples: no name overlap and sample counts differ.")
  }
}
write.table(alignment, file.path(out_dir, "sample_alignment.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

common_genes <- intersect(rownames(bulk), rownames(reference))
bulk <- bulk[common_genes, , drop = FALSE]
reference <- reference[common_genes, , drop = FALSE]
keep <- rowSums(bulk, na.rm = TRUE) > 0 & rowSums(reference, na.rm = TRUE) > 0
bulk <- bulk[keep, , drop = FALSE]
reference <- reference[keep, , drop = FALSE]
message("Matched genes used: ", nrow(bulk), "; cores: ", n_cores)

base_hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
base_hierarchy <- ImmuCellAI2.0:::validate_hierarchy(base_hierarchy, reference)
reference <- reference[, base_hierarchy$state_name, drop = FALSE]

valid_target_map <- target_map[names(target_map) %in% colnames(truth)]
valid_target_map <- valid_target_map[vapply(valid_target_map, function(x) any(x %in% colnames(reference)), logical(1))]
message("Comparable truth targets: ", length(valid_target_map))
write.table(
  data.frame(
    FlowCellType = names(valid_target_map),
    ReferenceStates = vapply(valid_target_map, paste, collapse = "+", FUN.VALUE = character(1))
  ),
  file.path(out_dir, "celltype_mapping_used.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

method_grid <- data.frame(
  Method = c(
    "ImmuCellAI2_tcell_VB",
    "ImmuCellAI2_tcell_VB_UNKNOWN",
    "ImmuCellAI2_flat_VB",
    "ImmuCellAI2_flat_VB_UNKNOWN"
  ),
  hierarchy.mode = c("tcell", "tcell", "flat", "flat"),
  inference.method = c("vb", "vb", "flat_vb", "flat_vb"),
  add.unknown = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

all_ih_metrics <- list()
all_ih_predictions <- list()

for (i in seq_len(nrow(method_grid))) {
  cfg <- method_grid[i, ]
  message("Running ", cfg$Method, " ...")
  hierarchy_i <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = cfg$hierarchy.mode)
  hierarchy_i <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_i, reference)
  ih_i <- deconvolve_bulk_matrix(
    bulk.mat = bulk,
    reference = reference,
    hierarchy = hierarchy_i,
    hierarchy.mode = cfg$hierarchy.mode,
    inference.method = cfg$inference.method,
    add.unknown = cfg$add.unknown,
    unknown.state = "UNKNOWN",
    unknown.mode = "mean",
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    n.cores = n_cores,
    seed = 123
  )
  state_names_i <- unique(c(hierarchy_i$state_name, "UNKNOWN"))
  pred_i <- ensure_samples_by_states(ih_i$state.fraction, colnames(bulk), state_names_i)
  all_ih_predictions[[cfg$Method]] <- pred_i
  all_ih_metrics[[cfg$Method]] <- evaluate_against_truth(pred_i, truth, valid_target_map, cfg$Method)
  write_deconvolution_outputs(ih_i, file.path(out_dir, cfg$Method))
  write.table(
    pred_i,
    file.path(out_dir, paste0(cfg$Method, "_state_fraction.txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  write.table(
    all_ih_metrics[[cfg$Method]],
    file.path(out_dir, paste0(cfg$Method, "_vs_truth_metrics.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message("Running BayesPrism chain.length=", chain_length, ", burn.in=", burn_in, "...")
reference_bp <- t(reference[, base_hierarchy$state_name, drop = FALSE])
mixture_bp <- t(bulk)
set.seed(123)
prism <- new.prism(
  reference = reference_bp,
  input.type = "count.matrix",
  cell.type.labels = base_hierarchy$major_lineage,
  cell.state.labels = base_hierarchy$state_name,
  key = NA_character_,
  mixture = mixture_bp,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
  pseudo.min = 1e-8
)
bp <- fig3_run_prism(
  prism,
  n.cores = n_cores,
  update.gibbs = FALSE,
  gibbs.control = list(
    chain.length = chain_length,
    burn.in = burn_in,
    thinning = 2,
    alpha = 1,
    seed = 123
  )
)
bp_pred <- get.fraction(bp, which.theta = "first", state.or.type = "state")
bp_pred <- ensure_samples_by_states(bp_pred, colnames(bulk), base_hierarchy$state_name)
write.table(bp_pred, file.path(out_dir, "bayesprism_state_fraction.txt"), sep = "\t", quote = FALSE, col.names = NA)

bp_metrics <- evaluate_against_truth(bp_pred, truth, valid_target_map, "BayesPrism_first_state_chain600_burn500")
write.table(bp_metrics, file.path(out_dir, "bayesprism_vs_truth_metrics.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

all_metrics <- do.call(rbind, c(all_ih_metrics, list(BayesPrism_first_state_chain600_burn500 = bp_metrics)))
write.table(all_metrics, file.path(out_dir, "all_methods_vs_truth_metrics.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

comparison_list <- lapply(names(all_ih_metrics), function(method_name) {
  x <- merge(
    all_ih_metrics[[method_name]][, c("CellType", "PredictedStates", "Pearson", "Spearman", "RMSE", "MAE", "Slope", "Intercept", "N_samples")],
    bp_metrics[, c("CellType", "Pearson", "Spearman", "RMSE", "MAE", "Slope", "Intercept", "N_samples")],
    by = "CellType",
    all = TRUE,
    suffixes = c("_ImmuCellAI2", "_BayesPrism")
  )
  x$ImmuCellAI2Method <- method_name
  x$DeltaPearson <- x$Pearson_ImmuCellAI2 - x$Pearson_BayesPrism
  x$DeltaSpearman <- x$Spearman_ImmuCellAI2 - x$Spearman_BayesPrism
  x$DeltaRMSE <- x$RMSE_ImmuCellAI2 - x$RMSE_BayesPrism
  x$DeltaMAE <- x$MAE_ImmuCellAI2 - x$MAE_BayesPrism
  x$BetterByPearson <- ifelse(
    x$DeltaPearson > 0,
    method_name,
    ifelse(x$DeltaPearson < 0, "BayesPrism", "Tie")
  )
  x
})
comparison_all <- do.call(rbind, comparison_list)
comparison_all <- comparison_all[order(comparison_all$ImmuCellAI2Method, -comparison_all$DeltaPearson), ]
write.table(
  comparison_all,
  file.path(out_dir, "per_celltype_correlation_difference_all_modes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_one <- function(metrics) {
  data.frame(
    N_targets = sum(!is.na(metrics$Pearson)),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE)
  )
}
summary_table <- do.call(rbind, lapply(split(all_metrics, all_metrics$Method), summary_one))
summary_table$Method <- rownames(summary_table)
rownames(summary_table) <- NULL
summary_table <- summary_table[, c("Method", setdiff(colnames(summary_table), "Method"))]
summary_table$GenesUsed <- nrow(bulk)
summary_table$SamplesUsed <- ncol(bulk)
write.table(summary_table, file.path(out_dir, "comparison_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

print(summary_table)
message("Done. Results written to: ", out_dir)
message("Main comparison table: ", file.path(out_dir, "per_celltype_correlation_difference_all_modes.txt"))
