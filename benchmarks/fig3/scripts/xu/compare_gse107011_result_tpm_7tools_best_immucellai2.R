# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
local_lib <- file.path(work_dir, "Rlib_test")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

stable_tmp <- file.path(work_dir, "Rtmp_immucellai2")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

bulk_file <- file.path(fig3_input_dir(), "result_tpm.csv")
truth_file <- file.path(fig3_input_dir(), "GSE107019_celltypeRatio.csv")
sample_id_file <- file.path(fig3_input_dir(), "GSE107019_sampleID.csv")
reference_file <- fig3_reference()
immune_gene_file <- fig3_markers()
out_dir <- file.path(fig3_work_dir(), "GSE107011_result_tpm_7tools_best_immucellai2")

n_cores <- fig3_cores()
chain_length <- 600
burn_in <- 500
cibersort_signature_genes <- 1000

# Cell types previously judged too broad or biologically ambiguous for this
# figure-style comparison. Edit this vector if you want to restore any target.
exclude_celltypes <- c(
  "B Ex", "Th17", "Tfh", "Th2", "T CD4 Memory", "B Naive", "MAIT",
  "Lymphoid", "Monocytes NC+I", "Tfh+Th", "Tregs", "Myeloid", "T CD4"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_pkgs <- c(
  "ImmuCellAI2.0", "BayesPrism", "CIBERSORT", "CITMIC", "ImmuCellAI",
  "MuSiC", "e1071", "limma"
)
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(ImmuCellAI2.0)
library(BayesPrism)
library(CIBERSORT)
library(CITMIC)
library(ImmuCellAI)
library(MuSiC)

patch_immucellai_gsva <- function() {
  ns <- asNamespace("ImmuCellAI")
  patched_sample_abundance <- function(sample, paper_marker, marker_exp, data_type, customer) {
    sam <- apply(sample, 2, as.numeric)
    row.names(sam) <- row.names(sample)
    tt <- intersect(row.names(sam), as.vector(unlist(paper_marker)))
    genes <- intersect(tt, row.names(marker_exp))
    sam_exp <- as.matrix(sam[genes, , drop = FALSE])
    colnames(sam_exp) <- colnames(sample)
    marker_exp <- marker_exp[genes, , drop = FALSE]
    marker_tag_mat <- c()
    for (cell in names(paper_marker)) {
      tag <- marker_tag(genes, as.vector(unlist(paper_marker[cell])))
      marker_tag_mat <- cbind(marker_tag_mat, tag)
    }
    row.names(marker_tag_mat) <- row.names(marker_exp)
    colnames(marker_tag_mat) <- names(paper_marker)
    exp_new <- apply(sam_exp, 2, function(x) sample_ratio(x, marker_exp, marker_tag_mat, data_type))
    if (requireNamespace("GSVA", quietly = TRUE) &&
        exists("ssgseaParam", envir = asNamespace("GSVA"))) {
      param <- GSVA::ssgseaParam(
        exprData = as.matrix(exp_new),
        geneSets = paper_marker,
        normalize = TRUE,
        verbose = FALSE
      )
      result <- GSVA::gsva(param)
    } else {
      result <- GSVA::gsva(as.matrix(exp_new), paper_marker, method = "ssgsea",
                           ssgsea.norm = TRUE, parallel.sz = 20)
    }
    if (ncol(result) < 3) {
      result[which(result < 0)] <- 0
    } else {
      result <- result - apply(result, 1, min)
    }
    compensation_matrix_num <- apply(compensation_matrix, 2, as.numeric)
    if (customer == 0) {
      row.names(compensation_matrix_num) <- row.names(compensation_matrix)
      result_norm <- compensation(result, compensation_matrix_num)
      if (ncol(result_norm) == 1) {
        InfiltrationScore <- sum(result_norm)
      } else if (nrow(result_norm) == 24) {
        InfiltrationScore <- apply(
          result_norm[c("Bcell", "CD4_T", "CD8_T", "DC", "Macrophage", "Monocyte", "Neutrophil", "NK"), ],
          2,
          sum
        )
      } else {
        InfiltrationScore <- apply(result_norm, 2, sum)
      }
      InfiltrationScore <- (InfiltrationScore / max(InfiltrationScore)) * 0.9
      result_mat <- rbind(result_norm, InfiltrationScore)
    } else {
      result_mat <- result
    }
    if (ncol(result_mat) > 1) {
      result_mat <- apply(result_mat, 1, function(x) round(x, 3))
    } else {
      result_mat <- t(round(result_mat, 3))
    }
    result_mat <- as.matrix(result_mat)
    if (ncol(result_mat) == 1) {
      colnames(result_mat) <- colnames(sample)
      T_FRE <- t(result_mat)
    } else {
      T_FRE <- result_mat
    }
    T_FRE
  }
  environment(patched_sample_abundance) <- ns
  was_locked <- bindingIsLocked("Sample_abundance_calculation", ns)
  if (was_locked) unlockBinding("Sample_abundance_calculation", ns)
  assign("Sample_abundance_calculation", patched_sample_abundance, envir = ns)
  if (was_locked) lockBinding("Sample_abundance_calculation", ns)
  invisible(TRUE)
}
patch_immucellai_gsva()

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
    file, sep = sep, header = FALSE, skip = 1, quote = "",
    comment.char = "", check.names = FALSE, fill = TRUE
  )
  genes <- as.character(dat[[1]])
  mat <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- if (length(header) == ncol(mat)) header else paste0("Sample", seq_len(ncol(mat)))
  rownames(mat) <- make.unique(genes)
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

read_expression_matrix <- function(file) {
  sep <- detect_sep(file)
  dat <- read.table(file, sep = sep, header = TRUE, row.names = 1, quote = "",
                    comment.char = "", check.names = FALSE)
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

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

read_gene_list <- function(file) {
  if (!file.exists(file)) return(character(0))
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[\t,; ]+"))
  x <- trimws(x)
  unique(x[nzchar(x)])
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

normalize_rows <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  rs <- rowSums(mat)
  rs[!is.finite(rs) | rs <= 0] <- eps
  sweep(mat, 1, rs, "/")
}

normalize_cols <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  cs <- colSums(mat)
  cs[!is.finite(cs) | cs <= 0] <- eps
  sweep(mat, 2, cs, "/")
}

normalize_nonnegative <- function(x) {
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  s <- sum(x)
  if (s <= 0) rep(1 / length(x), length(x)) else x / s
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

evaluate_against_truth <- function(pred, truth, target_map, method_name) {
  out <- lapply(names(target_map), function(cell) {
    if (!cell %in% colnames(truth)) return(NULL)
    states <- target_map[[cell]]
    pred_vec <- aggregate_prediction(pred, states)
    truth_vec <- truth[, cell]
    m <- metric_one(pred_vec, truth_vec)
    data.frame(
      Method = method_name,
      CellType = cell,
      PredictedStates = paste(intersect(states, colnames(pred)), collapse = "+"),
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

summary_eval <- function(metrics, method, genes.used, samples.used) {
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
    GenesUsed = genes.used,
    SamplesUsed = samples.used,
    ValidPearsonTargets = sum(!is.na(metrics$Pearson)),
    PenalizedMeanPearson = mean(ifelse(is.na(metrics$Pearson), -1, metrics$Pearson)),
    stringsAsFactors = FALSE
  )
}

state_target_map <- list(
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
  "mDCs" = "monoDC",
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

lm22_target_map <- list(
  "B Naive" = "B cells naive",
  "B NSM" = character(0),
  "B SM" = character(0),
  "Plasmablasts" = character(0),
  "T CD4 Naive" = "T cells CD4 naive",
  "Tregs" = "T cells regulatory (Tregs)",
  "Tfh" = "T cells follicular helper",
  "Th1" = character(0),
  "Th1/Th17" = character(0),
  "Th17" = character(0),
  "Th2" = character(0),
  "T CD8 Naive" = character(0),
  "T CD8 TE" = character(0),
  "MAIT" = character(0),
  "pDCs" = character(0),
  "mDCs" = c("Dendritic cells resting", "Dendritic cells activated"),
  "Monocytes C" = character(0),
  "Monocytes I" = character(0),
  "Neutrophils LD" = "Neutrophils",
  "Basophils LD" = character(0),
  "Tfh+Th1-17" = character(0),
  "T gd" = "T cells gamma delta",
  "Monocytes NC+I" = character(0),
  "Granulocytes LD" = character(0),
  "Tfh+Th" = c(
    "T cells CD4 memory resting", "T cells CD4 memory activated",
    "T cells follicular helper", "T cells regulatory (Tregs)"
  ),
  "T CD8 Memory" = character(0),
  "DCs" = c("Dendritic cells resting", "Dendritic cells activated"),
  "Monocytes C+I" = character(0),
  "T CD4 Memory" = c("T cells CD4 memory resting", "T cells CD4 memory activated"),
  "Monocytes" = "Monocytes",
  "T CD4" = c(
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)"
  ),
  "Myeloid Phagocytes" = c(
    "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated"
  ),
  "T Naive" = character(0),
  "T Memory" = character(0),
  "Lymphoid" = c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta",
    "NK cells resting", "NK cells activated"
  ),
  "Myeloid" = c(
    "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated",
    "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
  ),
  "Adaptive" = c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta"
  ),
  "Innate" = c(
    "NK cells resting", "NK cells activated", "Monocytes",
    "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated",
    "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
  )
)

citmic_target_map <- list(
  "B Naive" = "Naive B cell",
  "B NSM" = character(0),
  "B SM" = "Memory B cells",
  "Plasmablasts" = "Plasma cells",
  "T CD4 Naive" = "Naive CD4+ T cell",
  "Tregs" = c("Tregs", "nTreg", "iTregs"),
  "Tfh" = "Tfh",
  "Th1" = "Th1 cells",
  "Th1/Th17" = character(0),
  "Th17" = "Th17 cells",
  "Th2" = "Th2 cells",
  "T CD8 Naive" = "Naive CD8+ T cell",
  "T CD8 TE" = character(0),
  "MAIT" = "MAIT cells",
  "pDCs" = "pDCs",
  "mDCs" = "mDCs",
  "Monocytes C" = character(0),
  "Monocytes I" = character(0),
  "Neutrophils LD" = "Neutrophils",
  "Basophils LD" = "Basophils",
  "Tfh+Th1-17" = c("Tfh", "Th1 cells", "Th17 cells"),
  "T gd" = "Tgd cells",
  "Monocytes NC+I" = character(0),
  "Granulocytes LD" = c("Neutrophils", "Basophils", "Eosinophils"),
  "Tfh+Th" = c("T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells"),
  "T CD8 Memory" = c("CD8+ Tcm", "CD8+ Tem"),
  "DCs" = c("DCs", "mDCs", "pDCs"),
  "Monocytes C+I" = character(0),
  "T CD4 Memory" = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem"),
  "Monocytes" = "Monocytes",
  "T CD4" = c("CD4+ T cells", "T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells", "Tregs"),
  "Myeloid Phagocytes" = c("Monocytes", "Macrophages", "M1 Macrophages", "M2 Macrophages", "DCs", "mDCs", "pDCs"),
  "T Naive" = c("Naive CD4+ T cell", "Naive CD8+ T cell"),
  "T Memory" = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem", "CD8+ Tcm", "CD8+ Tem"),
  "Lymphoid" = c("B cells", "CD4+ T cells", "CD8+ T cells", "T helper cells", "NK cells", "NKT cells", "MAIT cells", "Tgd cells"),
  "Myeloid" = c("Monocytes", "Macrophages", "M1 Macrophages", "M2 Macrophages", "DCs", "mDCs", "pDCs", "Neutrophils", "Basophils", "Eosinophils", "Mast cell", "MDSCs"),
  "Adaptive" = c("B cells", "CD4+ T cells", "CD8+ T cells", "T helper cells", "Tregs"),
  "Innate" = c("MAIT cells", "Tgd cells", "NKT cells", "NK cells", "Monocytes", "Macrophages", "M1 Macrophages", "M2 Macrophages", "DCs", "mDCs", "pDCs", "Neutrophils", "Basophils", "Eosinophils", "Mast cell", "MDSCs")
)

immucellai_target_map <- list(
  "B Naive" = character(0),
  "B NSM" = character(0),
  "B SM" = character(0),
  "Plasmablasts" = character(0),
  "T CD4 Naive" = "CD4_naive",
  "Tregs" = c("nTreg", "iTreg", "Tr1"),
  "Tfh" = "Tfh",
  "Th1" = "Th1",
  "Th1/Th17" = character(0),
  "Th17" = "Th17",
  "Th2" = "Th2",
  "T CD8 Naive" = "CD8_naive",
  "T CD8 TE" = character(0),
  "MAIT" = "MAIT",
  "pDCs" = character(0),
  "mDCs" = character(0),
  "Monocytes C" = character(0),
  "Monocytes I" = character(0),
  "Neutrophils LD" = "Neutrophil",
  "Basophils LD" = character(0),
  "Tfh+Th1-17" = c("Tfh", "Th1", "Th17"),
  "T gd" = "Gamma_delta",
  "Monocytes NC+I" = character(0),
  "Granulocytes LD" = "Neutrophil",
  "Tfh+Th" = c("Tfh", "Th1", "Th17", "Th2"),
  "T CD8 Memory" = character(0),
  "DCs" = "DC",
  "Monocytes C+I" = character(0),
  "T CD4 Memory" = character(0),
  "Monocytes" = "Monocyte",
  "T CD4" = c("CD4_T", "CD4_naive", "Tfh", "Th1", "Th17", "Th2", "nTreg", "iTreg", "Tr1"),
  "Myeloid Phagocytes" = c("Monocyte", "Macrophage", "DC"),
  "T Naive" = c("CD4_naive", "CD8_naive"),
  "T Memory" = c("Central_memory", "Effector_memory"),
  "Lymphoid" = c("Bcell", "CD4_T", "CD8_T", "NK", "NKT", "MAIT", "Gamma_delta"),
  "Myeloid" = c("Monocyte", "Macrophage", "DC", "Neutrophil"),
  "Adaptive" = c("Bcell", "CD4_T", "CD8_T", "nTreg", "iTreg", "Tr1"),
  "Innate" = c("MAIT", "Gamma_delta", "NKT", "NK", "Monocyte", "Macrophage", "DC", "Neutrophil")
)

target_maps <- list(
  state = state_target_map,
  lm22 = lm22_target_map,
  citmic = citmic_target_map,
  immucellai = immucellai_target_map
)
target_maps <- lapply(target_maps, function(x) x[setdiff(names(x), exclude_celltypes)])

weighted_lm_nonnegative <- function(A, b, w = NULL) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  b[!is.finite(b)] <- 0
  if (is.null(w)) {
    Aw <- A
    bw <- b
  } else {
    w <- as.numeric(w)
    w[!is.finite(w) | w <= 0] <- median(w[w > 0 & is.finite(w)], na.rm = TRUE)
    w[!is.finite(w) | w <= 0] <- 1
    sw <- sqrt(w)
    Aw <- A * sw
    bw <- b * sw
  }
  coef <- tryCatch(as.numeric(stats::lm.fit(Aw, bw)$coefficients),
                   error = function(e) rep(1 / ncol(A), ncol(A)))
  normalize_nonnegative(coef)
}

dwls_one_sample <- function(A, b, n.iter = 6, eps = 1e-8) {
  b <- as.numeric(b)
  b[!is.finite(b)] <- 0
  if (sum(b) > 0) b <- b / sum(b)
  x <- weighted_lm_nonnegative(A, b)
  for (it in seq_len(n.iter)) {
    fitted <- as.numeric(A %*% x)
    resid <- abs(b - fitted)
    scale <- stats::quantile(resid, 0.75, na.rm = TRUE) + eps
    w <- 1 / (resid + scale)
    cap <- stats::quantile(w, 0.95, na.rm = TRUE)
    w <- pmin(w, cap)
    w <- w / mean(w, na.rm = TRUE)
    x.new <- weighted_lm_nonnegative(A, b, w)
    if (max(abs(x.new - x)) < 1e-7) {
      x <- x.new
      break
    }
    x <- x.new
  }
  normalize_nonnegative(x)
}

run_dwls_matrix <- function(bulk, reference) {
  ref.prob <- normalize_cols(reference)
  bulk.prob <- normalize_cols(bulk)
  run_one <- function(j) dwls_one_sample(ref.prob, bulk.prob[, j])
  fits <- if (.Platform$OS.type == "windows" && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("ref.prob", "bulk.prob", "dwls_one_sample",
                  "weighted_lm_nonnegative", "normalize_nonnegative"),
      envir = environment()
    )
    parallel::parLapply(cl, seq_len(ncol(bulk.prob)), function(j) dwls_one_sample(ref.prob, bulk.prob[, j]))
  } else {
    lapply(seq_len(ncol(bulk.prob)), run_one)
  }
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(bulk.prob)
  colnames(mat) <- colnames(ref.prob)
  mat
}

run_music_basic_matrix <- function(bulk, reference) {
  X <- normalize_cols(reference)
  Y <- normalize_cols(bulk)
  S <- rep(1, ncol(X))
  Sigma <- matrix(1e-8, nrow = nrow(X), ncol = ncol(X),
                  dimnames = list(rownames(X), colnames(X)))
  run_one <- function(j) {
    fit <- tryCatch(
      MuSiC::music.basic(Y = as.numeric(Y[, j]), X = X, S = S, Sigma = Sigma,
                         iter.max = 1000, nu = 1e-4, eps = 0.01),
      error = function(e) NULL
    )
    if (is.null(fit) || any(!is.finite(fit$p.weight))) weighted_lm_nonnegative(X, Y[, j]) else normalize_nonnegative(fit$p.weight)
  }
  fits <- lapply(seq_len(ncol(Y)), run_one)
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(Y)
  colnames(mat) <- colnames(X)
  mat
}

select_signature_genes <- function(reference, n = 1000) {
  ref <- log2(as.matrix(reference) + 1)
  state_var <- apply(ref, 1, stats::var, na.rm = TRUE)
  state_max <- apply(ref, 1, max, na.rm = TRUE)
  state_mean <- rowMeans(ref, na.rm = TRUE)
  score <- state_var * (state_max + 1e-8) / (state_mean + 1e-8)
  score[!is.finite(score)] <- 0
  genes <- names(sort(score, decreasing = TRUE))
  genes[seq_len(min(n, length(genes)))]
}

run_cibersort_default_lm22 <- function(bulk) {
  data(LM22, package = "CIBERSORT")
  genes <- intersect(rownames(LM22), rownames(bulk))
  fit <- CIBERSORT::cibersort(
    sig_matrix = LM22[genes, , drop = FALSE],
    mixture_file = bulk[genes, , drop = FALSE]
  )
  fit <- as.matrix(fit)
  fit[, intersect(colnames(LM22), colnames(fit)), drop = FALSE]
}

run_citmic_native_matrix <- function(bulk) {
  pred <- CITMIC::CITMIC(log2(as.data.frame(bulk) + 1), weighted = TRUE, base = 10,
                         damping = 0.90, cl.cores = n_cores)
  pred <- as.matrix(pred)
  if (!all(colnames(bulk) %in% rownames(pred)) && all(colnames(bulk) %in% colnames(pred))) pred <- t(pred)
  pred <- pred[colnames(bulk), , drop = FALSE]
  normalize_rows(pred)
}

run_immucellai_native_matrix <- function(bulk) {
  fit <- ImmuCellAI::ImmuCellAI_new(
    sample = as.data.frame(bulk),
    data_type = "RNA-seq",
    group_tag = FALSE,
    response_tag = FALSE,
    customer = 0
  )
  pred <- as.matrix(fit$Sample_abundance)
  pred <- pred[colnames(bulk), , drop = FALSE]
  normalize_rows(pred)
}

message("Reading input matrices...")
bulk <- collapse_duplicate_genes(read_bulk_result_tpm(bulk_file))
truth <- read_truth_matrix(truth_file)
reference <- collapse_duplicate_genes(read_expression_matrix(reference_file))

sample_map <- read_sample_id_map(sample_id_file)
if (!is.null(sample_map)) {
  sample_map <- sample_map[sample_map$BulkSample %in% colnames(bulk) &
                             sample_map$TruthSample %in% rownames(truth), , drop = FALSE]
  if (nrow(sample_map) < 3) stop("Sample ID map has fewer than 3 overlapping samples.")
  sample_map <- sample_map[match(intersect(colnames(bulk), sample_map$BulkSample), sample_map$BulkSample), , drop = FALSE]
  bulk <- bulk[, sample_map$BulkSample, drop = FALSE]
  truth <- truth[sample_map$TruthSample, , drop = FALSE]
  rownames(truth) <- sample_map$BulkSample
  alignment <- sample_map
} else {
  common_samples <- intersect(colnames(bulk), rownames(truth))
  if (length(common_samples) >= 3) {
    bulk <- bulk[, common_samples, drop = FALSE]
    truth <- truth[common_samples, , drop = FALSE]
    alignment <- data.frame(BulkSample = common_samples, TruthSample = common_samples)
  } else if (ncol(bulk) == nrow(truth)) {
    alignment <- data.frame(BulkSample = colnames(bulk), TruthSample = rownames(truth))
    rownames(truth) <- colnames(bulk)
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

base_hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
base_hierarchy <- ImmuCellAI2.0:::validate_hierarchy(base_hierarchy, reference)
reference <- reference[, base_hierarchy$state_name, drop = FALSE]

immune_genes <- intersect(read_gene_list(immune_gene_file), rownames(reference))
if (length(immune_genes) >= 100) {
  bulk_ih <- bulk[immune_genes, , drop = FALSE]
  reference_ih <- reference[immune_genes, , drop = FALSE]
} else {
  bulk_ih <- bulk
  reference_ih <- reference
}

valid_truth_cells <- colnames(truth)
target_maps <- lapply(target_maps, function(x) x[names(x) %in% valid_truth_cells])
write.table(
  data.frame(
    CellType = names(target_maps$state),
    ReferenceStates = vapply(target_maps$state, paste, collapse = "+", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "celltype_mapping_state_tools.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

run_status <- data.frame(Method = character(), Status = character(), Message = character())
add_status <- function(method, status, msg = "") {
  run_status <<- rbind(run_status, data.frame(Method = method, Status = status, Message = msg))
  write.table(run_status, file.path(out_dir, "run_status.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

method_grid <- data.frame(
  Method = c("ImmuCellAI2_tcell_VB", "ImmuCellAI2_tcell_VB_UNKNOWN",
             "ImmuCellAI2_flat_VB", "ImmuCellAI2_flat_VB_UNKNOWN"),
  hierarchy.mode = c("tcell", "tcell", "flat", "flat"),
  inference.method = c("vb", "vb", "flat_vb", "flat_vb"),
  add.unknown = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

ih_metrics <- list()
ih_summary <- list()
for (i in seq_len(nrow(method_grid))) {
  cfg <- method_grid[i, ]
  method <- cfg$Method
  fraction_file <- file.path(out_dir, paste0(method, "_state_fraction.txt"))
  add_status(method, "RUNNING", "")
  if (file.exists(fraction_file)) {
    pred <- read.delim(fraction_file, row.names = 1, check.names = FALSE)
    pred <- ensure_samples_by_states(pred, colnames(bulk))
    add_status(method, "REUSED", fraction_file)
  } else {
    hierarchy_i <- create_default_53_hierarchy(colnames(reference_ih), hierarchy.mode = cfg$hierarchy.mode)
    hierarchy_i <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_i, reference_ih)
    fit <- deconvolve_bulk_matrix(
      bulk.mat = bulk_ih,
      reference = reference_ih,
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
    pred <- ensure_samples_by_states(fit$state.fraction, colnames(bulk))
    write_deconvolution_outputs(fit, file.path(out_dir, method))
    write.table(pred, fraction_file, sep = "\t", quote = FALSE, col.names = NA)
    add_status(method, "OK", "")
  }
  ih_metrics[[method]] <- evaluate_against_truth(pred, truth, target_maps$state, method)
  ih_summary[[method]] <- summary_eval(ih_metrics[[method]], method, nrow(bulk_ih), ncol(bulk))
  write.table(ih_metrics[[method]], file.path(out_dir, paste0(method, "_vs_truth_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

immucellai2_mode_summary <- do.call(rbind, ih_summary)
immucellai2_mode_summary <- immucellai2_mode_summary[order(immucellai2_mode_summary$PenalizedMeanPearson, decreasing = TRUE), ]
rownames(immucellai2_mode_summary) <- NULL
write.table(immucellai2_mode_summary, file.path(out_dir, "immucellai2_mode_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

best_ih_method <- immucellai2_mode_summary$Method[1]
best_ih_metrics <- ih_metrics[[best_ih_method]]
best_ih_metrics$Method <- paste0("ImmuCellAI2_best_", best_ih_method)
writeLines(best_ih_method, file.path(out_dir, "best_immucellai2_mode.txt"))

method_metrics <- list()
method_metrics[[paste0("ImmuCellAI2_best_", best_ih_method)]] <- best_ih_metrics

add_or_reuse <- function(method, fraction_file, run_expr, target_map, genes.used = nrow(bulk)) {
  add_status(method, "RUNNING", "")
  if (file.exists(fraction_file)) {
    pred <- read.delim(fraction_file, row.names = 1, check.names = FALSE)
    pred <- ensure_samples_by_states(pred, colnames(bulk))
    add_status(method, "REUSED", fraction_file)
  } else {
    pred <- run_expr()
    pred <- ensure_samples_by_states(pred, colnames(bulk))
    write.table(pred, fraction_file, sep = "\t", quote = FALSE, col.names = NA)
    add_status(method, "OK", "")
  }
  metrics <- evaluate_against_truth(pred, truth, target_map, method)
  write.table(metrics, file.path(out_dir, paste0(method, "_vs_truth_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  method_metrics[[method]] <<- metrics
  invisible(metrics)
}

add_or_reuse(
  "BayesPrism_first_state_chain600_burn500",
  file.path(out_dir, "BayesPrism_first_state_chain600_burn500_state_fraction.txt"),
  function() {
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
      gibbs.control = list(chain.length = chain_length, burn.in = burn_in,
                           thinning = 2, alpha = 1, seed = 123)
    )
    get.fraction(bp, which.theta = "first", state.or.type = "state")
  },
  target_maps$state
)

add_or_reuse("DWLS_weighted_lm", file.path(out_dir, "DWLS_weighted_lm_fraction.txt"),
             function() run_dwls_matrix(bulk, reference), target_maps$state)
add_or_reuse("MuSiC_basic", file.path(out_dir, "MuSiC_basic_fraction.txt"),
             function() run_music_basic_matrix(bulk, reference), target_maps$state)
add_or_reuse("CIBERSORT_default_LM22", file.path(out_dir, "CIBERSORT_default_LM22_fraction.txt"),
             function() run_cibersort_default_lm22(bulk), target_maps$lm22, genes.used = 547)
add_or_reuse("CITMIC_native", file.path(out_dir, "CITMIC_native_fraction.txt"),
             function() run_citmic_native_matrix(bulk), target_maps$citmic)
add_or_reuse("ImmuCellAI_native", file.path(out_dir, "ImmuCellAI_native_fraction.txt"),
             function() run_immucellai_native_matrix(bulk), target_maps$immucellai)

all_metrics <- do.call(rbind, method_metrics)
write.table(all_metrics, file.path(out_dir, "all_methods_per_celltype_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- do.call(rbind, lapply(names(method_metrics), function(method) {
  summary_eval(method_metrics[[method]], method, nrow(bulk), ncol(bulk))
}))
summary_table <- summary_table[order(summary_table$PenalizedMeanPearson, decreasing = TRUE), ]
rownames(summary_table) <- NULL
write.table(summary_table, file.path(out_dir, "comparison_summary_7tools.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pearson_wide <- reshape(
  all_metrics[, c("CellType", "Method", "Pearson")],
  idvar = "CellType",
  timevar = "Method",
  direction = "wide"
)
colnames(pearson_wide) <- sub("^Pearson\\.", "", colnames(pearson_wide))
write.table(pearson_wide, file.path(out_dir, "per_celltype_pearson_7tools_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

spearman_wide <- reshape(
  all_metrics[, c("CellType", "Method", "Spearman")],
  idvar = "CellType",
  timevar = "Method",
  direction = "wide"
)
colnames(spearman_wide) <- sub("^Spearman\\.", "", colnames(spearman_wide))
write.table(spearman_wide, file.path(out_dir, "per_celltype_spearman_7tools_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(immucellai2_mode_summary)
print(summary_table)
message("Best ImmuCellAI2 mode: ", best_ih_method)
message("Done. Results written to: ", out_dir)
