# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
out_dir <- file.path(work_dir, "GSE107011_result_tpm_7tools_best_immucellai2")
truth_file <- file.path(fig3_input_dir(), "GSE107019_celltypeRatio.csv")
sample_id_file <- file.path(fig3_input_dir(), "GSE107019_sampleID.csv")

read_truth_matrix <- function(file) {
  dat <- read.csv(file, skip = 1, header = TRUE, check.names = FALSE,
                  quote = "\"", comment.char = "")
  sample_col <- grep("^Sample", colnames(dat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sample_col)) sample_col <- colnames(dat)[1]
  sample_names <- as.character(dat[[sample_col]])
  dat <- dat[, setdiff(colnames(dat), sample_col), drop = FALSE]
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- sample_names
  mat[is.na(mat)] <- NA_real_
  if (max(mat, na.rm = TRUE) > 1.5) mat <- mat / 100
  mat
}

read_sample_id_map <- function(file) {
  if (!file.exists(file)) return(NULL)
  dat <- read.csv(file, header = FALSE, check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  if (ncol(dat) < 2) stop("sample_id_file must contain at least two columns.")
  out <- data.frame(
    BulkSample = trimws(as.character(dat[[1]])),
    TruthSample = trimws(as.character(dat[[2]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$BulkSample) & nzchar(out$TruthSample), , drop = FALSE]
  out[!duplicated(out$BulkSample), , drop = FALSE]
}

ensure_samples_by_states <- function(pred, sample_names) {
  pred <- as.matrix(pred)
  if (all(sample_names %in% rownames(pred))) {
    pred <- pred[sample_names, , drop = FALSE]
  } else if (all(sample_names %in% colnames(pred))) {
    pred <- t(pred[, sample_names, drop = FALSE])
  }
  pred
}

aggregate_by_cols <- function(mat, cols) {
  if (length(cols) == 0) return(rep(NA_real_, nrow(mat)))
  cols <- intersect(cols, colnames(mat))
  if (length(cols) == 0) return(rep(NA_real_, nrow(mat)))
  rowSums(mat[, cols, drop = FALSE], na.rm = TRUE)
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
    return(c(Pearson = NA, Spearman = NA, RMSE = NA, MAE = NA,
             Slope = NA, Intercept = NA, N = 0))
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

summary_eval <- function(metrics, method, samples.used) {
  data.frame(
    Method = method,
    N_targets = nrow(metrics),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MedianSpearman = median(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MedianRMSE = median(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE),
    MedianMAE = median(metrics$MAE, na.rm = TRUE),
    SamplesUsed = samples.used,
    ValidPearsonTargets = sum(!is.na(metrics$Pearson)),
    PenalizedMeanPearson = mean(ifelse(is.na(metrics$Pearson), -1, metrics$Pearson)),
    stringsAsFactors = FALSE
  )
}

targets <- data.frame(
  FullCellType = c(
    "memory B cell", "Basophils", "B naive", "T CD4 Naive",
    "CD8Tcm", "T CD8 Naive", "pDCs", "mDCs", "Plasmablasts",
    "classical Monocytes", "intermediate Monocytes", "Th1", "Tfh",
    "Th1/Th17", "Neutrophils", "gdT", "Treg", "MAIT", "Monocytes",
    "CD4 T", "CD8 T", "DC", "Myeloid Phagocytes", "Adaptive"
  ),
  CellType = c(
    "MBC", "Basophil", "Bn", "CD4Tn",
    "CD8Tcm", "CD8Tn", "pDC", "mDC", "PB",
    "cMo", "intMo", "Th1", "Tfh", "Th1.17",
    "Neutrophil", "γδT", "Treg", "MAIT", "Monocyte",
    "CD4T", "CD8T", "DC", "MyeloidPhagocytes", "Adaptive"
  ),
  stringsAsFactors = FALSE
)

truth_map <- list(
  "MBC" = c("B NSM", "B SM"),
  "Basophil" = "Basophils LD",
  "Bex" = "B Ex",
  "Bn" = "B Naive",
  "CD4Tn" = "T CD4 Naive",
  "CD8Tcm" = "T CD8 CM",
  "CD8Tn" = "T CD8 Naive",
  "pDC" = "pDCs",
  "mDC" = "mDCs",
  "PB" = "Plasmablasts",
  "cMo" = "Monocytes C",
  "intMo" = "Monocytes I",
  "Th1" = "Th1",
  "Tfh" = "Tfh",
  "Th1.17" = "Th1/Th17",
  "Neutrophil" = "Neutrophils LD",
  "γδT" = "T gd",
  "Treg" = "Tregs",
  "MAIT" = "MAIT",
  "Monocyte" = "Monocytes",
  "CD4T" = "T CD4",
  "CD8T" = "T CD8",
  "DC" = "DCs",
  "MyeloidPhagocytes" = "Myeloid Phagocytes",
  "Adaptive" = "Adaptive"
)

state_map <- list(
  "MBC" = c("FOB", "MBC"),
  "Basophil" = "Basophil",
  "Bex" = "Bex",
  "Bn" = "Bnaive",
  "CD4Tn" = "CD4Tn",
  "CD8Tcm" = "CD8Tcm",
  "CD8Tn" = "CD8Tn",
  "pDC" = "pDC",
  "mDC" = "monoDC",
  "PB" = "PB",
  "cMo" = "cMo",
  "intMo" = "intMo",
  "Th1" = "Th1",
  "Tfh" = "Tfh",
  "Th1.17" = "Th1/17",
  "Neutrophil" = "Neutrophil",
  "γδT" = "gdT",
  "Treg" = "Treg",
  "MAIT" = "MAIT",
  "Monocyte" = c("cMo", "intMo", "ncMo"),
  "CD4T" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
             "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  "CD8T" = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex"),
  "DC" = c("cDC1", "cDC2", "pDC", "monoDC"),
  "MyeloidPhagocytes" = c("cMo", "intMo", "ncMo", "M0", "M1", "M2", "TAM",
                          "cDC1", "cDC2", "pDC", "monoDC"),
  "Adaptive" = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC",
                 "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
                 "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm",
                 "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg")
)

lm22_map <- list(
  "MBC" = "B cells memory",
  "Basophil" = character(0),
  "Bex" = character(0),
  "Bn" = "B cells naive",
  "CD4Tn" = "T cells CD4 naive",
  "CD8Tcm" = character(0),
  "CD8Tn" = character(0),
  "pDC" = character(0),
  "mDC" = c("Dendritic cells resting", "Dendritic cells activated"),
  "PB" = character(0),
  "cMo" = character(0),
  "intMo" = character(0),
  "Th1" = character(0),
  "Tfh" = "T cells follicular helper",
  "Th1.17" = character(0),
  "Neutrophil" = "Neutrophils",
  "纬未T" = "T cells gamma delta",
  "Treg" = "T cells regulatory (Tregs)",
  "MAIT" = character(0),
  "Monocyte" = "Monocytes",
  "CD4T" = c(
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)"
  ),
  "CD8T" = "T cells CD8",
  "DC" = c("Dendritic cells resting", "Dendritic cells activated"),
  "MyeloidPhagocytes" = c(
    "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated"
  ),
  "Adaptive" = c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta"
  )
)

citmic_map <- list(
  "MBC" = c("Memory B cells", "Classswitched memory B cells"),
  "Basophil" = "Basophils",
  "Bex" = character(0),
  "Bn" = "Naive B cell",
  "CD4Tn" = "Naive CD4+ T cell",
  "CD8Tcm" = "CD8+ Tcm",
  "CD8Tn" = "Naive CD8+ T cell",
  "pDC" = "pDCs",
  "mDC" = "mDCs",
  "PB" = "Plasma cells",
  "cMo" = character(0),
  "intMo" = character(0),
  "Th1" = "Th1 cells",
  "Tfh" = "Tfh",
  "Th1.17" = character(0),
  "Neutrophil" = "Neutrophils",
  "γδT" = "Tgd cells",
  "Treg" = c("Tregs", "nTreg", "iTregs"),
  "MAIT" = "MAIT cells",
  "Monocyte" = "Monocytes",
  "CD4T" = "CD4+ T cells",
  "CD8T" = "CD8+ T cells",
  "DC" = "DCs",
  "MyeloidPhagocytes" = c("Monocytes", "Macrophages", "M1 Macrophages", "M2 Macrophages",
                          "DCs", "mDCs", "pDCs"),
  "Adaptive" = c("B cells", "CD4+ T cells", "CD8+ T cells", "Tregs")
)

immucellai_map <- list(
  "MBC" = character(0),
  "Basophil" = character(0),
  "Bex" = character(0),
  "Bn" = character(0),
  "CD4Tn" = "CD4_naive",
  "CD8Tcm" = character(0),
  "CD8Tn" = "CD8_naive",
  "pDC" = character(0),
  "mDC" = character(0),
  "PB" = character(0),
  "cMo" = character(0),
  "intMo" = character(0),
  "Th1" = "Th1",
  "Tfh" = "Tfh",
  "Th1.17" = character(0),
  "Neutrophil" = "Neutrophil",
  "γδT" = "Gamma_delta",
  "Treg" = c("nTreg", "iTreg", "Tr1"),
  "MAIT" = "MAIT",
  "Monocyte" = "Monocyte",
  "CD4T" = c("CD4_T", "CD4_naive", "Th1", "Th17", "Th2", "Tfh", "nTreg", "iTreg", "Tr1"),
  "CD8T" = c("CD8_T", "CD8_naive", "Cytotoxic", "Exhausted"),
  "DC" = "DC",
  "MyeloidPhagocytes" = c("Monocyte", "Macrophage", "DC"),
  "Adaptive" = c("Bcell", "CD4_T", "CD8_T", "nTreg", "iTreg", "Tr1")
)

evaluate_against_truth <- function(pred, truth, pred_map, method_name) {
  out <- lapply(seq_len(nrow(targets)), function(i) {
    short <- targets$CellType[i]
    truth_cols <- truth_map[[short]]
    pred_cols <- pred_map[[short]]
    pred_vec <- aggregate_by_cols(pred, pred_cols)
    truth_vec <- aggregate_by_cols(truth, truth_cols)
    m <- metric_one(pred_vec, truth_vec)
    data.frame(
      Method = method_name,
      FullCellType = targets$FullCellType[i],
      CellType = short,
      TruthColumns = paste(intersect(truth_cols, colnames(truth)), collapse = "+"),
      PredictedStates = paste(intersect(pred_cols, colnames(pred)), collapse = "+"),
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
  do.call(rbind, out)
}

message("Reading truth and aligning samples...")
truth <- read_truth_matrix(truth_file)
sample_map <- read_sample_id_map(sample_id_file)
if (!is.null(sample_map)) {
  sample_map <- sample_map[sample_map$TruthSample %in% rownames(truth), , drop = FALSE]
  truth <- truth[sample_map$TruthSample, , drop = FALSE]
  rownames(truth) <- sample_map$BulkSample
}

read_fraction <- function(file) {
  x <- read.delim(file, row.names = 1, check.names = FALSE)
  ensure_samples_by_states(x, rownames(truth))
}

method_files <- c(
  "ImmuCellAI2_tcell_VB" = "ImmuCellAI2_tcell_VB_state_fraction.txt",
  "ImmuCellAI2_tcell_VB_UNKNOWN" = "ImmuCellAI2_tcell_VB_UNKNOWN_state_fraction.txt",
  "ImmuCellAI2_flat_VB" = "ImmuCellAI2_flat_VB_state_fraction.txt",
  "ImmuCellAI2_flat_VB_UNKNOWN" = "ImmuCellAI2_flat_VB_UNKNOWN_state_fraction.txt",
  "BayesPrism_first_state_chain600_burn500" = "BayesPrism_first_state_chain600_burn500_state_fraction.txt",
  "DWLS_weighted_lm" = "DWLS_weighted_lm_fraction.txt",
  "MuSiC_basic" = "MuSiC_basic_fraction.txt",
  "CIBERSORT_default_LM22" = "CIBERSORT_default_LM22_fraction.txt",
  "CITMIC_native" = "CITMIC_native_fraction.txt",
  "ImmuCellAI_native" = "ImmuCellAI_native_fraction.txt"
)

preds <- lapply(method_files, function(f) read_fraction(file.path(out_dir, f)))
common_samples <- Reduce(intersect, c(list(rownames(truth)), lapply(preds, rownames)))
if (length(common_samples) < 3) stop("Fewer than 3 common samples across truth and prediction files.")
truth <- truth[common_samples, , drop = FALSE]
preds <- lapply(preds, function(x) x[common_samples, , drop = FALSE])

message("Recalculating ImmuCellAI2 mode metrics...")
ih_methods <- names(method_files)[grepl("^ImmuCellAI2_", names(method_files))]
ih_metrics <- lapply(ih_methods, function(method) {
  evaluate_against_truth(preds[[method]], truth, state_map, method)
})
names(ih_metrics) <- ih_methods
ih_summary <- do.call(rbind, lapply(ih_methods, function(method) {
  summary_eval(ih_metrics[[method]], method, nrow(truth))
}))
ih_summary <- ih_summary[order(ih_summary$PenalizedMeanPearson, decreasing = TRUE), ]
rownames(ih_summary) <- NULL
best_ih_method <- ih_summary$Method[1]
writeLines(best_ih_method, file.path(out_dir, "best_immucellai2_mode.txt"))
write.table(ih_summary, file.path(out_dir, "immucellai2_mode_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (method in ih_methods) {
  write.table(ih_metrics[[method]], file.path(out_dir, paste0(method, "_vs_truth_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Recalculating seven-tool metrics...")
method_metrics <- list()
best_metrics <- ih_metrics[[best_ih_method]]
best_metrics$Method <- paste0("ImmuCellAI2_best_", best_ih_method)
method_metrics[[best_metrics$Method[1]]] <- best_metrics
method_metrics[["BayesPrism_first_state_chain600_burn500"]] <-
  evaluate_against_truth(preds[["BayesPrism_first_state_chain600_burn500"]], truth, state_map,
                         "BayesPrism_first_state_chain600_burn500")
method_metrics[["DWLS_weighted_lm"]] <-
  evaluate_against_truth(preds[["DWLS_weighted_lm"]], truth, state_map, "DWLS_weighted_lm")
method_metrics[["MuSiC_basic"]] <-
  evaluate_against_truth(preds[["MuSiC_basic"]], truth, state_map, "MuSiC_basic")
method_metrics[["CIBERSORT_default_LM22"]] <-
  evaluate_against_truth(preds[["CIBERSORT_default_LM22"]], truth, lm22_map,
                         "CIBERSORT_default_LM22")
method_metrics[["CITMIC_native"]] <-
  evaluate_against_truth(preds[["CITMIC_native"]], truth, citmic_map, "CITMIC_native")
method_metrics[["ImmuCellAI_native"]] <-
  evaluate_against_truth(preds[["ImmuCellAI_native"]], truth, immucellai_map, "ImmuCellAI_native")

for (method in setdiff(names(method_metrics), best_metrics$Method[1])) {
  write.table(method_metrics[[method]], file.path(out_dir, paste0(method, "_vs_truth_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

all_metrics <- do.call(rbind, method_metrics)
write.table(all_metrics, file.path(out_dir, "all_methods_per_celltype_metricsGSE107011.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_metrics, file.path(out_dir, "all_methods_per_celltype_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- do.call(rbind, lapply(names(method_metrics), function(method) {
  summary_eval(method_metrics[[method]], method, nrow(truth))
}))
summary_table <- summary_table[order(summary_table$PenalizedMeanPearson, decreasing = TRUE), ]
rownames(summary_table) <- NULL
write.table(summary_table, file.path(out_dir, "comparison_summary_7tools.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pearson_wide <- reshape(
  all_metrics[, c("FullCellType", "CellType", "Method", "Pearson")],
  idvar = c("FullCellType", "CellType"),
  timevar = "Method",
  direction = "wide"
)
colnames(pearson_wide) <- sub("^Pearson\\.", "", colnames(pearson_wide))
write.table(pearson_wide, file.path(out_dir, "per_celltype_pearson_7tools_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

spearman_wide <- reshape(
  all_metrics[, c("FullCellType", "CellType", "Method", "Spearman")],
  idvar = c("FullCellType", "CellType"),
  timevar = "Method",
  direction = "wide"
)
colnames(spearman_wide) <- sub("^Spearman\\.", "", colnames(spearman_wide))
write.table(spearman_wide, file.path(out_dir, "per_celltype_spearman_7tools_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  data.frame(
    FullCellType = targets$FullCellType,
    CellType = targets$CellType,
    TruthColumns = vapply(targets$CellType, function(x) paste(intersect(truth_map[[x]], colnames(truth)), collapse = "+"), character(1)),
    ReferenceStates = vapply(targets$CellType, function(x) paste(state_map[[x]], collapse = "+"), character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "celltype_mapping_state_tools.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

run_status <- data.frame(
  Method = names(method_files),
  Status = "REUSED_FRACTION_RECALCULATED_METRICS",
  Message = unname(method_files),
  stringsAsFactors = FALSE
)
write.table(run_status, file.path(out_dir, "run_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(ih_summary)
print(summary_table)
message("Best ImmuCellAI2 mode: ", best_ih_method)
message("Done. Recalculated requested GSE107011 cell types in: ", out_dir)
