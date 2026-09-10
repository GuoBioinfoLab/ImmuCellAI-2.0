# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
out_dir <- file.path(work_dir, "GSE107011_adjusted_7tools_flatvb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

truth_file <- file.path(fig3_input_dir(), "GSE107019_celltypeRatio.csv")
sample_id_file <- file.path(fig3_input_dir(), "GSE107019_sampleID.csv")
fraction_dir <- file.path(work_dir, "GSE107011_result_tpm_7tools_best_immucellai2")
old_immucellai2_bayesprism_dir <- file.path(fig3_input_dir(), "GSE107019_compare_ImmuCellAI2_BayesPrism")

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
  CellType = c(
    "T CD8 Memory", "T Memory", "Monocytes C", "Monocytes I",
    "T CD4 Naive", "Basophils LD", "T CD8 TE", "Innate",
    "Granulocytes LD", "Monocytes", "memory B cell", "Adaptive",
    "T gd", "DCs", "Plasmablasts",
    "T helper", "Neutrophils LD", "T CD8 Naive",
    "T CD4", "mDC", "pDCs"
  ),
  TruthKey = c(
    "T CD8 Memory", "T Memory", "Monocytes C", "Monocytes I",
    "T CD4 Naive", "Basophils LD", "T CD8 TE", "Innate",
    "Granulocytes LD", "Monocytes", "memory B cell", "Adaptive",
    "T gd", "DCs", "Plasmablasts",
    "Tfh+Th", "Neutrophils LD", "T CD8 Naive",
    "T CD4", "mDC", "pDCs"
  ),
  stringsAsFactors = FALSE
)

truth_map <- list(
  "T CD8 Memory" = "T CD8 Memory",
  "T Memory" = "T Memory",
  "Monocytes C" = "Monocytes C",
  "Th17" = "Th17",
  "Monocytes I" = "Monocytes I",
  "T CD4 Naive" = "T CD4 Naive",
  "Th1" = "Th1",
  "Basophils LD" = "Basophils LD",
  "T CD8 TE" = "T CD8 TE",
  "Innate" = "Innate",
  "MAIT" = "MAIT",
  "Granulocytes LD" = "Granulocytes LD",
  "Monocytes" = "Monocytes",
  "memory B cell" = c("B NSM", "B SM"),
  "Adaptive" = "Adaptive",
  "T gd" = "T gd",
  "Th1/Th17" = "Th1/Th17",
  "DCs" = "DCs",
  "Plasmablasts" = "Plasmablasts",
  "T CD4 Memory" = "T CD4 Memory",
  "Tfh+Th" = "Tfh+Th",
  "Neutrophils LD" = "Neutrophils LD",
  "Tfh+Th1-17" = "Tfh+Th1-17",
  "T CD8 Naive" = "T CD8 Naive",
  "T CD4" = "T CD4",
  "mDC" = "mDCs",
  "pDCs" = "pDCs"
)

state_map <- list(
  "T CD8 Memory" = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
  "T Memory" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
  "Monocytes C" = "cMo",
  "Th17" = "Th17",
  "Monocytes I" = "intMo",
  "T CD4 Naive" = "CD4Tn",
  "Th1" = "Th1",
  "Basophils LD" = "Basophil",
  "T CD8 TE" = "CD8Temra",
  "Innate" = c("MAIT", "gdT", "NKT", "cNK", "NKreg", "ILC1", "ILC2", "ILC3", "cMo", "intMo", "ncMo", "M0", "M1", "M2", "TAM", "cDC1", "cDC2", "pDC", "monoDC", "Neutrophil", "Basophil", "Eosinophil", "Mast cell", "MDSC"),
  "MAIT" = "MAIT",
  "Granulocytes LD" = c("Neutrophil", "Basophil", "Eosinophil"),
  "Monocytes" = c("cMo", "intMo", "ncMo"),
  "memory B cell" = c("FOB", "MBC"),
  "Adaptive" = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC", "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  "T gd" = "gdT",
  "Th1/Th17" = "Th1/17",
  "DCs" = c("cDC1", "cDC2", "pDC", "monoDC"),
  "Plasmablasts" = "PB",
  "T CD4 Memory" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
  "Tfh+Th" = c("Tfh", "Th1", "Th1/17", "Th17", "Th2"),
  "Neutrophils LD" = "Neutrophil",
  "Tfh+Th1-17" = c("Tfh", "Th1", "Th1/17", "Th17"),
  "T CD8 Naive" = "CD8Tn",
  "T CD4" = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  "mDC" = "monoDC",
  "pDCs" = "pDC"
)

lm22_map <- list(
  "T CD8 Memory" = character(0),
  "T Memory" = character(0),
  "Monocytes C" = character(0),
  "Th17" = character(0),
  "Monocytes I" = character(0),
  "T CD4 Naive" = "T cells CD4 naive",
  "Th1" = character(0),
  "Basophils LD" = character(0),
  "T CD8 TE" = character(0),
  "Innate" = c(
    "T cells gamma delta", "NK cells resting", "NK cells activated",
    "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated",
    "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
  ),
  "MAIT" = character(0),
  "Granulocytes LD" = character(0),
  "Monocytes" = "Monocytes",
  "memory B cell" = "B cells memory",
  "Adaptive" = c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta"
  ),
  "T gd" = "T cells gamma delta",
  "Th1/Th17" = character(0),
  "DCs" = c("Dendritic cells resting", "Dendritic cells activated"),
  "Plasmablasts" = character(0),
  "T CD4 Memory" = c("T cells CD4 memory resting", "T cells CD4 memory activated"),
  "Tfh+Th" = c(
    "T cells CD4 memory resting", "T cells CD4 memory activated",
    "T cells follicular helper", "T cells regulatory (Tregs)"
  ),
  "Neutrophils LD" = "Neutrophils",
  "Tfh+Th1-17" = character(0),
  "T CD8 Naive" = character(0),
  "T CD4" = c(
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)"
  ),
  "mDC" = c("Dendritic cells resting", "Dendritic cells activated"),
  "pDCs" = character(0)
)

citmic_map <- list(
  "T CD8 Memory" = c("CD8+ Tcm", "CD8+ Tem"),
  "T Memory" = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem", "CD8+ Tcm", "CD8+ Tem"),
  "Monocytes C" = character(0),
  "Th17" = "Th17 cells",
  "Monocytes I" = character(0),
  "T CD4 Naive" = "Naive CD4+ T cell",
  "Th1" = "Th1 cells",
  "Basophils LD" = "Basophils",
  "T CD8 TE" = character(0),
  "Innate" = c("MAIT cells", "Tgd cells", "NKT cells", "NK cells", "Monocytes", "Macrophages", "M1 Macrophages", "M2 Macrophages", "DCs", "mDCs", "pDCs", "Neutrophils", "Basophils", "Eosinophils", "Mast cell", "MDSCs"),
  "MAIT" = "MAIT cells",
  "Granulocytes LD" = c("Neutrophils", "Basophils", "Eosinophils"),
  "Monocytes" = "Monocytes",
  "memory B cell" = c("Memory B cells", "Classswitched memory B cells"),
  "Adaptive" = c("B cells", "CD4+ T cells", "CD8+ T cells", "T helper cells", "Tregs"),
  "T gd" = "Tgd cells",
  "Th1/Th17" = character(0),
  "DCs" = c("DCs", "mDCs", "pDCs"),
  "Plasmablasts" = "Plasma cells",
  "T CD4 Memory" = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem"),
  "Tfh+Th" = c("T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells"),
  "Neutrophils LD" = "Neutrophils",
  "Tfh+Th1-17" = c("Tfh", "Th1 cells", "Th17 cells"),
  "T CD8 Naive" = "Naive CD8+ T cell",
  "T CD4" = c("CD4+ T cells", "T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells", "Tregs"),
  "mDC" = "mDCs",
  "pDCs" = "pDCs"
)

immucellai_map <- list(
  "T CD8 Memory" = character(0),
  "T Memory" = c("Central_memory", "Effector_memory"),
  "Monocytes C" = character(0),
  "Th17" = "Th17",
  "Monocytes I" = character(0),
  "T CD4 Naive" = "CD4_naive",
  "Th1" = "Th1",
  "Basophils LD" = character(0),
  "T CD8 TE" = character(0),
  "Innate" = c("MAIT", "Gamma_delta", "NKT", "NK", "Monocyte", "Macrophage", "DC", "Neutrophil"),
  "MAIT" = "MAIT",
  "Granulocytes LD" = "Neutrophil",
  "Monocytes" = "Monocyte",
  "memory B cell" = character(0),
  "Adaptive" = c("Bcell", "CD4_T", "CD8_T", "nTreg", "iTreg", "Tr1"),
  "T gd" = "Gamma_delta",
  "Th1/Th17" = character(0),
  "DCs" = "DC",
  "Plasmablasts" = character(0),
  "T CD4 Memory" = character(0),
  "Tfh+Th" = c("Tfh", "Th1", "Th17", "Th2"),
  "Neutrophils LD" = "Neutrophil",
  "Tfh+Th1-17" = c("Tfh", "Th1", "Th17"),
  "T CD8 Naive" = "CD8_naive",
  "T CD4" = c("CD4_T", "CD4_naive", "Tfh", "Th1", "Th17", "Th2", "nTreg", "iTreg", "Tr1"),
  "mDC" = character(0),
  "pDCs" = character(0)
)

evaluate_against_truth <- function(pred, truth, pred_map, method_name) {
  out <- lapply(seq_len(nrow(targets)), function(i) {
    label <- targets$CellType[i]
    key <- targets$TruthKey[i]
    truth_cols <- truth_map[[key]]
    pred_cols <- pred_map[[key]]
    pred_vec <- aggregate_by_cols(pred, pred_cols)
    truth_vec <- aggregate_by_cols(truth, truth_cols)
    m <- metric_one(pred_vec, truth_vec)
    data.frame(
      Method = method_name,
      CellType = label,
      TruthKey = key,
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

truth <- read_truth_matrix(truth_file)
sample_map <- read_sample_id_map(sample_id_file)
sample_map <- sample_map[sample_map$TruthSample %in% rownames(truth), , drop = FALSE]
truth <- truth[sample_map$TruthSample, , drop = FALSE]
rownames(truth) <- sample_map$BulkSample

method_files <- c(
  "ImmuCellAI2_flat_VB" = file.path(old_immucellai2_bayesprism_dir, "ImmuCellAI2_flat_VB_state_fraction.txt"),
  "BayesPrism_first_state_chain600_burn500" = file.path(old_immucellai2_bayesprism_dir, "bayesprism_state_fraction.txt"),
  "DWLS_weighted_lm" = file.path(fraction_dir, "DWLS_weighted_lm_fraction.txt"),
  "MuSiC_basic" = file.path(fraction_dir, "MuSiC_basic_fraction.txt"),
  "CIBERSORT_default_LM22" = file.path(fraction_dir, "CIBERSORT_default_LM22_fraction.txt"),
  "CITMIC_native" = file.path(fraction_dir, "CITMIC_native_fraction.txt"),
  "ImmuCellAI_native" = file.path(fraction_dir, "ImmuCellAI_native_fraction.txt")
)

read_fraction <- function(file) {
  x <- read.delim(file, row.names = 1, check.names = FALSE)
  ensure_samples_by_states(x, rownames(truth))
}

preds <- lapply(method_files, read_fraction)
common_samples <- Reduce(intersect, c(list(rownames(truth)), lapply(preds, rownames)))
truth <- truth[common_samples, , drop = FALSE]
preds <- lapply(preds, function(x) x[common_samples, , drop = FALSE])

method_metrics <- list(
  "ImmuCellAI2_flat_VB" = evaluate_against_truth(preds[["ImmuCellAI2_flat_VB"]], truth, state_map, "ImmuCellAI2_flat_VB"),
  "BayesPrism_first_state_chain600_burn500" = evaluate_against_truth(preds[["BayesPrism_first_state_chain600_burn500"]], truth, state_map, "BayesPrism_first_state_chain600_burn500"),
  "DWLS_weighted_lm" = evaluate_against_truth(preds[["DWLS_weighted_lm"]], truth, state_map, "DWLS_weighted_lm"),
  "MuSiC_basic" = evaluate_against_truth(preds[["MuSiC_basic"]], truth, state_map, "MuSiC_basic"),
  "CIBERSORT_default_LM22" = evaluate_against_truth(preds[["CIBERSORT_default_LM22"]], truth, lm22_map, "CIBERSORT_default_LM22"),
  "CITMIC_native" = evaluate_against_truth(preds[["CITMIC_native"]], truth, citmic_map, "CITMIC_native"),
  "ImmuCellAI_native" = evaluate_against_truth(preds[["ImmuCellAI_native"]], truth, immucellai_map, "ImmuCellAI_native")
)

all_metrics <- do.call(rbind, method_metrics)
write.table(all_metrics, file.path(out_dir, "GSE107011_adjusted_7tools_per_celltype_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- do.call(rbind, lapply(names(method_metrics), function(method) {
  summary_eval(method_metrics[[method]], method, nrow(truth))
}))
summary_table <- summary_table[order(summary_table$PenalizedMeanPearson, decreasing = TRUE), ]
rownames(summary_table) <- NULL
write.table(summary_table, file.path(out_dir, "GSE107011_adjusted_7tools_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pearson_wide <- reshape(
  all_metrics[, c("CellType", "Method", "Pearson")],
  idvar = "CellType",
  timevar = "Method",
  direction = "wide"
)
colnames(pearson_wide) <- sub("^Pearson\\.", "", colnames(pearson_wide))
write.table(pearson_wide, file.path(out_dir, "GSE107011_adjusted_7tools_pearson_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package RColorBrewer is required.")
if (!requireNamespace("magrittr", quietly = TRUE)) stop("Package magrittr is required.")
if (!requireNamespace("stringr", quietly = TRUE)) stop("Package stringr is required.")
library(ggplot2)
library(magrittr)
if (requireNamespace("scRNAtoolVis", quietly = TRUE)) library(scRNAtoolVis)

native_geom_jjpie <- exists("geom_jjpie", mode = "function")
method_labels <- c(
  ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  DWLS_weighted_lm = "DWLS",
  MuSiC_basic = "MuSiC",
  CITMIC_native = "CITMIC",
  CIBERSORT_default_LM22 = "CIBERSORT",
  ImmuCellAI_native = "ImmuCellAI"
)
method_order <- c("ImmuCellAI 2.0", "BayesPrism", "DWLS", "MuSiC", "CITMIC", "CIBERSORT", "ImmuCellAI")

display_label_map <- c(
  "T CD8 Memory" = "CD8Tmemory",
  "T Memory" = "Tmemory",
  "Monocytes C" = "CMonocyte",
  "Monocytes I" = "IMonocyte",
  "T CD4 Naive" = "CD4Tn",
  "Basophils LD" = "Basophil",
  "T CD8 TE" = "CD8Temra",
  "Granulocytes LD" = "Granulocyte",
  "memory B cell" = "Bmemory",
  "T gd" = "gdT",
  "Plasmablasts" = "plasma",
  "T CD4 Memory" = "CD4Tmemory",
  "T helper" = "T helper",
  "T CD8 Naive" = "CD8Tn"
)

plot_data <- all_metrics
plot_data$Dataset <- "GSE107011"
plot_data$method <- unname(method_labels[plot_data$Method])
plot_data$cellType <- ifelse(
  plot_data$CellType %in% names(display_label_map),
  unname(display_label_map[plot_data$CellType]),
  plot_data$CellType
)
plot_data$correlation <- suppressWarnings(as.numeric(plot_data$Pearson))
plot_data$status <- "valid"
plot_data$status[is.na(plot_data$correlation) & trimws(plot_data$PredictedStates) == ""] <- "not_available"
plot_data$status[is.na(plot_data$correlation) & trimws(plot_data$PredictedStates) != ""] <- "constant_prediction"

write.table(plot_data, file.path(out_dir, "GSE107011_adjusted_7tools_pie_plot_data.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

theme_blue <- theme(
  plot.title = element_text(size = 13, face = "bold", color = "darkred", hjust = 0, lineheight = 1.2),
  plot.subtitle = element_text(size = 13, face = "bold", color = "grey30", lineheight = 1.2, hjust = 0),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "gray80", linewidth = 1, linetype = "dashed"),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(vjust = 1, face = "bold", size = 14, color = "darkred"),
  axis.title.y = element_text(size = 14, face = "bold", color = "darkred"),
  axis.text.x = element_text(size = 12, colour = "black"),
  legend.title = element_text(size = 12, colour = "black"),
  axis.text.y = element_text(size = 12, colour = "black"),
  legend.text = element_text(size = 12, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  legend.key = element_blank(),
  strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 1),
  strip.text = element_text(size = 18, colour = "black")
)

make_circle_poly <- function(cx, cy, r, n = 96) {
  a <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}
make_wedge_poly <- function(cx, cy, r, frac, start = pi / 2, n = 96) {
  frac <- min(max(abs(frac), 0), 1)
  if (!is.finite(frac) || frac <= 0) return(NULL)
  a <- seq(start, start - 2 * pi * frac, length.out = max(3, ceiling(n * frac) + 1))
  data.frame(x = c(cx, cx + r * cos(a), cx), y = c(cy, cy + r * sin(a), cy))
}
if (!native_geom_jjpie) {
  geom_jjpie <- function(data, mapping = NULL, width = 0.65, ...) {
    if (is.null(data) || nrow(data) == 0) return(list())
    df <- data
    r <- min(width / 2, 0.48)
    bg_list <- vector("list", nrow(df))
    wedge_list <- list()
    for (i in seq_len(nrow(df))) {
      bg <- make_circle_poly(df$x_id[i], df$y_id[i], r = r)
      bg$idx <- i
      bg_list[[i]] <- bg
      if (is.finite(df$correlation[i])) {
        wg <- make_wedge_poly(df$x_id[i], df$y_id[i], r = r, frac = min(abs(df$correlation[i]), 1))
        if (!is.null(wg)) {
          wg$idx <- i
          wg$correlation <- df$correlation[i]
          wedge_list[[length(wedge_list) + 1]] <- wg
        }
      }
    }
    bg_poly <- do.call(rbind, bg_list)
    wedge_poly <- if (length(wedge_list) > 0) do.call(rbind, wedge_list) else data.frame()
    list(
      geom_polygon(data = bg_poly, aes(x = x, y = y, group = idx), inherit.aes = FALSE,
                   fill = "white", color = "#1f1f1f", linewidth = 0.65),
      geom_polygon(data = wedge_poly, aes(x = x, y = y, group = idx, fill = correlation), inherit.aes = FALSE,
                   color = "#1f1f1f", linewidth = 0.55)
    )
  }
}

cell_order <- sort(unique(as.character(plot_data$cellType)))
plot_data$method <- factor(plot_data$method, levels = rev(method_order))
plot_data$cellType <- factor(plot_data$cellType, levels = cell_order)
plot_data$x_id <- as.numeric(plot_data$cellType)
plot_data$y_id <- as.numeric(plot_data$method)
constant_df <- plot_data[plot_data$status == "constant_prediction", , drop = FALSE]
not_available_df <- plot_data[plot_data$status == "not_available", , drop = FALSE]

if (native_geom_jjpie) {
  p <- ggplot(plot_data, aes(x = cellType, y = method, fill = correlation)) +
    geom_point(data = subset(plot_data, !is.na(correlation)), aes(x = cellType, y = method),
               inherit.aes = FALSE, shape = 21, size = 2.2, fill = "white",
               color = "#1f1f1f", stroke = 1.05) +
    geom_jjpie(data = subset(plot_data, !is.na(correlation)), aes(piefill = correlation),
               width = 0.65, color = "#1f1f1f", linewidth = 0.55) +
    geom_point(data = constant_df, aes(x = cellType, y = method),
               inherit.aes = FALSE, shape = 4, size = 1.5, color = "black", stroke = 1.8) +
    geom_text(data = not_available_df, aes(x = cellType, y = method),
              inherit.aes = FALSE, label = "/", color = "black", size = 2.8, fontface = "bold") +
    facet_grid(. ~ Dataset, scales = "free", space = "free") +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(), limits = c(-1, 1), name = "Pearson") +
    scale_x_discrete(expand = expansion(add = c(0.55, 0.55))) +
    scale_y_discrete(limits = unique(plot_data$method) %>% stringr::str_sort() %>% rev(),
                     expand = expansion(add = c(0.55, 0.55))) +
    labs(x = NULL, y = NULL, caption = "Undefined glyphs: x = constant prediction; slash = not available") +
    theme_blue +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"))
} else {
  axis_x <- data.frame(x = seq_along(cell_order), label = cell_order)
  axis_y <- data.frame(y = seq_along(levels(plot_data$method)), label = levels(plot_data$method))
  slash_df <- data.frame()
  if (nrow(not_available_df) > 0) {
    slash_df <- do.call(rbind, lapply(seq_len(nrow(not_available_df)), function(i) {
      data.frame(x = not_available_df$x_id[i] - 0.13, y = not_available_df$y_id[i] - 0.13,
                 xend = not_available_df$x_id[i] + 0.13, yend = not_available_df$y_id[i] + 0.13, idx = i)
    }))
  }
  legend_status <- data.frame(
    x = c(Inf, Inf), y = c(Inf, Inf),
    status = factor(c("constant prediction", "not available"), levels = c("constant prediction", "not available"))
  )
  p <- ggplot(plot_data, aes(x = x_id, y = y_id, fill = correlation)) +
    geom_jjpie(data = subset(plot_data, !is.na(correlation)), aes(piefill = correlation), width = 0.65) +
    geom_point(data = constant_df, aes(x = x_id, y = y_id), inherit.aes = FALSE,
               shape = 4, size = 1.5, color = "black", stroke = 1.8) +
    geom_segment(data = slash_df, aes(x = x, y = y, xend = xend, yend = yend, group = idx),
                 inherit.aes = FALSE, color = "black", linewidth = 0.95, lineend = "round") +
    geom_point(data = legend_status, aes(x = x, y = y, shape = status), inherit.aes = FALSE, alpha = 0) +
    facet_grid(. ~ Dataset, scales = "free", space = "free") +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral")), limits = c(-1, 1), name = "Pearson") +
    scale_shape_manual(values = c(`constant prediction` = "x", `not available` = "/"), name = "Undefined") +
    guides(shape = guide_legend(override.aes = list(alpha = 1, color = "black", size = 3))) +
    scale_x_continuous(breaks = axis_x$x, labels = axis_x$label, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = axis_y$y, labels = axis_y$label, expand = expansion(mult = c(0.08, 0.08))) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = NULL, caption = "Undefined glyphs: x = constant prediction; slash = not available") +
    theme_blue +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"))
}

width_high <- length(cell_order) * 0.88 + 5.2
height_high <- length(method_order) * 0.84 + 3.0
ggsave(file.path(out_dir, "GSE107011_adjusted_7tools_pie_high.pdf"), p, width = width_high, height = height_high, limitsize = FALSE)
ggsave(file.path(out_dir, "GSE107011_adjusted_7tools_pie_low.pdf"), p, width = width_high, height = height_high, limitsize = FALSE)
ggsave(file.path(out_dir, "GSE107011_adjusted_7tools_pie_high.png"), p, width = width_high, height = height_high, dpi = 220, limitsize = FALSE)

print(summary_table)
message("Saved output directory: ", out_dir)
