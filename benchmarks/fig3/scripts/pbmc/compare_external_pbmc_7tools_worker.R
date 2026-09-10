# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

# Worker script: compare seven deconvolution methods against external PBMC truth.
#
# Methods:
#   1. ImmuCellAI2_tcell_VB
#   2. BayesPrism_first_state_chain600_burn500
#   3. DWLS_native_ref_top1000
#   4. MuSiC_native_music.basic
#   5. CITMIC_native
#   6. CIBERSORT_default_LM22
#   7. ImmuCellAI_native
#
# Run:

options(stringsAsFactors = FALSE)

workspace_dir <- fig3_work_dir()
local_lib <- file.path(workspace_dir, "Rlib_test")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

stable_tmp <- file.path(workspace_dir, "Rtmp_immucellai2")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

if (!exists("dataset_name", inherits = FALSE)) stop("dataset_name must be supplied by the runner.")
if (!exists("bulk_file", inherits = FALSE)) stop("bulk_file must be supplied by the runner.")
if (!exists("truth_file", inherits = FALSE)) stop("truth_file must be supplied by the runner.")
if (!exists("reference_file", inherits = FALSE)) {
  reference_file <- fig3_reference()
}
if (!exists("immune_gene_file", inherits = FALSE)) {
  immune_gene_file <- fig3_markers()
}
if (!exists("out_dir", inherits = FALSE)) stop("out_dir must be supplied by the runner.")
if (!exists("immucellai2_mode", inherits = FALSE)) immucellai2_mode <- "flat"
if (!exists("immucellai2_add_unknown", inherits = FALSE)) immucellai2_add_unknown <- FALSE
if (!exists("immucellai2_method_name", inherits = FALSE)) {
  immucellai2_method_name <- paste0(
    "ImmuCellAI2_", immucellai2_mode, "_VB",
    if (immucellai2_add_unknown) "_UNKNOWN" else ""
  )
}

if (!exists("sample_limit", inherits = FALSE)) sample_limit <- NA
if (!exists("n_cores", inherits = FALSE)) n_cores <- min(8, parallel::detectCores(logical = TRUE))

if (!exists("run_immucellai2", inherits = FALSE)) run_immucellai2 <- TRUE
if (!exists("run_bayesprism", inherits = FALSE)) run_bayesprism <- TRUE
if (!exists("run_dwls", inherits = FALSE)) run_dwls <- TRUE
if (!exists("run_music", inherits = FALSE)) run_music <- TRUE
if (!exists("run_citmic", inherits = FALSE)) run_citmic <- TRUE
if (!exists("run_cibersort", inherits = FALSE)) run_cibersort <- TRUE
if (!exists("run_immucellai", inherits = FALSE)) run_immucellai <- TRUE
if (!exists("cibersort_signature_genes", inherits = FALSE)) cibersort_signature_genes <- 1000

if (!exists("chain_length", inherits = FALSE)) chain_length <- 600
if (!exists("burn_in", inherits = FALSE)) burn_in <- 500

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_pkgs <- c(
  "ImmuCellAI2.0", "BayesPrism", "CIBERSORT", "CITMIC", "ImmuCellAI",
  "MuSiC", "nnls", "e1071"
)
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(ImmuCellAI2.0)
if (isTRUE(run_bayesprism)) library(BayesPrism)
if (isTRUE(run_citmic)) library(CITMIC)
if (isTRUE(run_cibersort)) library(CIBERSORT)
if (isTRUE(run_immucellai)) library(ImmuCellAI)
if (isTRUE(run_music)) library(MuSiC)

patch_immucellai_gsva <- function() {
  if (!isTRUE(run_immucellai)) return(invisible(FALSE))
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
      } else {
        if (nrow(result_norm) == 24) {
          InfiltrationScore <- apply(
            result_norm[c("Bcell", "CD4_T", "CD8_T", "DC", "Macrophage", "Monocyte", "Neutrophil", "NK"), ],
            2,
            sum
          )
        } else {
          InfiltrationScore <- apply(result_norm, 2, sum)
        }
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

read_gene_list <- function(file) {
  if (!file.exists(file)) return(character(0))
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[\t,; ]+"))
  x <- trimws(x)
  unique(x[nzchar(x)])
}

normalize_cols <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  cs <- colSums(mat)
  cs[!is.finite(cs) | cs <= 0] <- eps
  sweep(mat, 2, cs, "/")
}

normalize_rows <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  rs <- rowSums(mat)
  rs[!is.finite(rs) | rs <= 0] <- eps
  sweep(mat, 1, rs, "/")
}

normalize_nonnegative <- function(x) {
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  s <- sum(x)
  if (s <= 0) rep(1 / length(x), length(x)) else x / s
}

run_dwls_matrix <- function(bulk, reference) {
  genes <- select_signature_genes(reference, cibersort_signature_genes)
  S <- normalize_cols(reference[genes, , drop = FALSE])
  B <- normalize_cols(bulk[genes, , drop = FALSE])
  dwls.cores <- max(1, min(as.integer(n_cores), 8))
  run_one <- function(j) {
    as.numeric(DWLS::solveDampenedWLS(S = S, B = B[, j]))
  }
  if (dwls.cores <= 1) {
    fits <- lapply(seq_len(ncol(B)), run_one)
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(dwls.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(DWLS))
    parallel::clusterExport(
      cl,
      varlist = c("S", "B"),
      envir = environment()
    )
    fits <- parallel::parLapply(cl, seq_len(ncol(B)), function(j) {
      as.numeric(DWLS::solveDampenedWLS(S = S, B = B[, j]))
    })
  } else {
    fits <- parallel::mclapply(seq_len(ncol(B)), run_one, mc.cores = dwls.cores)
  }
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(B)
  colnames(mat) <- colnames(S)
  mat
}

run_music_basic_matrix <- function(bulk, reference) {
  X <- normalize_cols(reference)
  Y <- normalize_cols(bulk)
  S <- rep(1, ncol(X))
  Sigma <- matrix(1e-8, nrow = nrow(X), ncol = ncol(X),
                  dimnames = list(rownames(X), colnames(X)))
  music.cores <- max(1, min(as.integer(n_cores), 8))
  run_one <- function(j) {
    fit <- tryCatch(
      MuSiC::music.basic(Y = as.numeric(Y[, j]), X = X, S = S, Sigma = Sigma,
                         iter.max = 1000, nu = 1e-4, eps = 0.01),
      error = function(e) NULL
    )
    if (is.null(fit) || any(!is.finite(fit$p.weight))) stop("MuSiC::music.basic failed.")
    normalize_nonnegative(fit$p.weight)
  }
  if (music.cores <= 1) {
    fits <- lapply(seq_len(ncol(Y)), run_one)
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(music.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(MuSiC))
    parallel::clusterExport(
      cl,
      varlist = c("X", "Y", "S", "Sigma", "normalize_nonnegative"),
      envir = environment()
    )
    fits <- parallel::parLapply(cl, seq_len(ncol(Y)), function(j) {
      fit <- tryCatch(
        MuSiC::music.basic(Y = as.numeric(Y[, j]), X = X, S = S, Sigma = Sigma,
                           iter.max = 1000, nu = 1e-4, eps = 0.01),
        error = function(e) NULL
      )
      if (is.null(fit) || any(!is.finite(fit$p.weight))) stop("MuSiC::music.basic failed.")
      normalize_nonnegative(fit$p.weight)
    })
  } else {
    fits <- parallel::mclapply(seq_len(ncol(Y)), run_one, mc.cores = music.cores)
  }
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

run_cibersort_lm22_matrix <- function(bulk) {
  data(LM22, package = "CIBERSORT")
  genes <- intersect(rownames(LM22), rownames(bulk))
  fit <- CIBERSORT::cibersort(
    sig_matrix = as.matrix(LM22[genes, , drop = FALSE]),
    mixture_file = as.matrix(bulk[genes, , drop = FALSE])
  )
  fit <- as.matrix(fit)
  fraction_cols <- intersect(colnames(LM22), colnames(fit))
  fit[, fraction_cols, drop = FALSE]
}

run_citmic_native_matrix <- function(bulk) {
  citmic.cores <- max(1, min(as.integer(n_cores), 8))
  citmic.input <- log2(as.matrix(bulk) + 1)
  pred <- CITMIC::CITMIC(as.data.frame(citmic.input), weighted = TRUE, base = 10,
                         damping = 0.90, cl.cores = citmic.cores)
  pred <- as.matrix(pred)
  if (!all(colnames(bulk) %in% rownames(pred)) && all(colnames(bulk) %in% colnames(pred))) {
    pred <- t(pred)
  }
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

extract_variable_sample_info <- function(samples) {
  parsed <- regexec("^(.+)_([0-9]+)$", samples)
  pieces <- regmatches(samples, parsed)
  ok <- lengths(pieces) == 3
  if (!all(ok)) {
    stop("Cannot parse sample names: ", paste(utils::head(samples[!ok], 10), collapse = ", "))
  }
  cell_type <- vapply(pieces, `[`, character(1), 2)
  index <- as.integer(vapply(pieces, `[`, character(1), 3))
  data.frame(
    Sample = samples,
    CellType = cell_type,
    StepIndex = index,
    TrueFraction = index * 0.005,
    stringsAsFactors = FALSE
  )
}

state_target_mapper <- function(target) {
  mapper <- list(
    Bcells = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC"),
    CD4Tcm = "CD4Tcm",
    CD4Tem = "CD4Tem",
    CD4Temra = "CD4Temra",
    CD4Tfh = "Tfh",
    CD4Tnaive = "CD4Tn",
    CD4Trm = "CD4Trm",
    CD8Tcm = "CD8Tcm",
    CD8Tem = "CD8Tem",
    CD8Temra = "CD8Temra",
    CD8Tnaive = "CD8Tn",
    CD8Trm = "CD8Trm",
    MAIT = "MAIT",
    Monocyte = c("cMo", "intMo", "ncMo"),
    NK = c("cNK", "NKreg"),
    NKT = "NKT",
    Th1 = "Th1",
    Th17 = "Th17",
    monoDC = "monoDC",
    `γδT` = "gdT",
    CD4Tmemory = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
    CD4Tn = "CD4Tn",
    CD8Tmemory = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
    CD8Tn = "CD8Tn",
    CMonocyte = "cMo",
    FOB = "FOB",
    GC_B = "BGC",
    ILC3 = "ILC3",
    NMonocyte = "ncMo",
    TAM = "TAM",
    Thelper = c("Tfh", "Th1", "Th1/17", "Th17", "Th2"),
    Treg = "Treg",
    cDC1 = "cDC1",
    cDC2 = "cDC2",
    cytotoxicNK = "cNK",
    exhausted_T = "Tex",
    macrophage = c("M0", "M1", "M2"),
    mast_cell = "Mast cell",
    pDC = "pDC",
    plasma = c("PB", "PC"),
    regulatoryNK = "NKreg"
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

citmic_target_mapper <- function(target) {
  mapper <- list(
    Bcells = "B cells",
    CD4Tcm = "CD4+ Tcm",
    CD4Tem = "CD4+ Tem",
    CD4Temra = character(0),
    CD4Tfh = "Tfh",
    CD4Tnaive = "Naive CD4+ T cell",
    CD4Trm = character(0),
    CD8Tcm = "CD8+ Tcm",
    CD8Tem = "CD8+ Tem",
    CD8Temra = character(0),
    CD8Tnaive = "Naive CD8+ T cell",
    CD8Trm = character(0),
    MAIT = "MAIT cells",
    Monocyte = "Monocytes",
    NK = "NK cells",
    NKT = "NKT cells",
    Th1 = "Th1 cells",
    Th17 = "Th17 cells",
    monoDC = character(0),
    `γδT` = "Tgd cells",
    CD4Tmemory = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem"),
    CD4Tn = "Naive CD4+ T cell",
    CD8Tmemory = c("CD8+ Tcm", "CD8+ Tem"),
    CD8Tn = "Naive CD8+ T cell",
    CMonocyte = character(0),
    FOB = character(0),
    GC_B = character(0),
    ILC3 = character(0),
    NMonocyte = character(0),
    TAM = character(0),
    Thelper = c("T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells"),
    Treg = c("Tregs", "nTreg", "iTregs"),
    cDC1 = character(0),
    cDC2 = character(0),
    cytotoxicNK = "NK CD56dim cells",
    exhausted_T = "Exhausted T cells",
    macrophage = c("Macrophages", "M1 Macrophages", "M2 Macrophages"),
    mast_cell = "Mast cell",
    pDC = "pDCs",
    plasma = "Plasma cells",
    regulatoryNK = "NK CD56bright cells"
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

immucellai_target_mapper <- function(target) {
  mapper <- list(
    Bcells = "Bcell",
    CD4Tcm = character(0),
    CD4Tem = character(0),
    CD4Temra = character(0),
    CD4Tfh = "Tfh",
    CD4Tnaive = "CD4_naive",
    CD4Trm = character(0),
    CD8Tcm = character(0),
    CD8Tem = character(0),
    CD8Temra = character(0),
    CD8Tnaive = "CD8_naive",
    CD8Trm = character(0),
    MAIT = "MAIT",
    Monocyte = "Monocyte",
    NK = "NK",
    NKT = "NKT",
    Th1 = "Th1",
    Th17 = "Th17",
    monoDC = character(0),
    `γδT` = "Gamma_delta",
    CD4Tmemory = character(0),
    CD4Tn = "CD4_naive",
    CD8Tmemory = character(0),
    CD8Tn = "CD8_naive",
    CMonocyte = character(0),
    FOB = character(0),
    GC_B = character(0),
    ILC3 = character(0),
    NMonocyte = character(0),
    TAM = character(0),
    Thelper = c("Th1", "Th2", "Th17", "Tfh"),
    Treg = c("nTreg", "iTreg", "Tr1"),
    cDC1 = character(0),
    cDC2 = character(0),
    cytotoxicNK = character(0),
    exhausted_T = "Exhausted",
    macrophage = "Macrophage",
    mast_cell = character(0),
    pDC = character(0),
    plasma = character(0),
    regulatoryNK = character(0)
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

evaluate_prediction <- function(pred.mat, sample.info, target.mapper) {
  pred.mat <- as.matrix(pred.mat)
  if (!all(sample.info$Sample %in% rownames(pred.mat)) &&
      all(sample.info$Sample %in% colnames(pred.mat))) {
    pred.mat <- t(pred.mat)
  }
  pred.mat <- pred.mat[sample.info$Sample, , drop = FALSE]
  targets <- unique(sample.info$CellType)
  rows <- lapply(targets, function(ct) {
    states <- intersect(target.mapper(ct)[[ct]], colnames(pred.mat))
    idx <- sample.info$CellType == ct
    truth <- sample.info$TrueFraction[idx]
    pred <- if (length(states) == 0) rep(NA_real_, sum(idx)) else rowSums(pred.mat[idx, states, drop = FALSE], na.rm = TRUE)
    ok <- is.finite(pred) & is.finite(truth)
    fit <- if (sum(ok) >= 3 && stats::sd(pred[ok]) > 0 && stats::sd(truth[ok]) > 0) {
      stats::coef(stats::lm(pred[ok] ~ truth[ok]))
    } else {
      c(`(Intercept)` = NA_real_, truth = NA_real_)
    }
    data.frame(
      CellType = ct,
      PredictedStates = paste(states, collapse = "+"),
      Pearson = safe_cor(pred, truth, "pearson"),
      Spearman = safe_cor(pred, truth, "spearman"),
      RMSE = if (sum(ok) > 0) sqrt(mean((pred[ok] - truth[ok])^2)) else NA_real_,
      MAE = if (sum(ok) > 0) mean(abs(pred[ok] - truth[ok])) else NA_real_,
      Slope = unname(fit[2]),
      Intercept = unname(fit[1]),
      N = sum(ok),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_eval <- function(metrics, method, genes.used) {
  data.frame(
    Method = method,
    Level = "external_ground_truth",
    N_targets = nrow(metrics),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = stats::median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MedianSpearman = stats::median(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MedianRMSE = stats::median(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE),
    MedianMAE = stats::median(metrics$MAE, na.rm = TRUE),
    MeanSlope = mean(metrics$Slope, na.rm = TRUE),
    MedianSlope = stats::median(metrics$Slope, na.rm = TRUE),
    GenesUsed = genes.used,
    ValidPearsonTargets = sum(!is.na(metrics$Pearson)),
    PenalizedMeanPearson = mean(ifelse(is.na(metrics$Pearson), -1, metrics$Pearson)),
    stringsAsFactors = FALSE
  )
}

external_state_map <- list(
  Neutrophils = "Neutrophil",
  Monocytes = c("cMo", "intMo", "ncMo"),
  CD8T = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex"),
  CD4T = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
           "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  B = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC"),
  NK = c("cNK", "NKreg")
)
external_citmic_map <- list(
  Neutrophils = "Neutrophils",
  Monocytes = "Monocytes",
  CD8T = "CD8+ T cells",
  CD4T = "CD4+ T cells",
  B = "B cells",
  NK = "NK cells"
)
external_immucellai_map <- list(
  Neutrophils = "Neutrophil",
  Monocytes = "Monocyte",
  CD8T = "CD8_T",
  CD4T = "CD4_T",
  B = "Bcell",
  NK = "NK"
)
external_lm22_map <- list(
  Neutrophils = "Neutrophils",
  Monocytes = "Monocytes",
  CD8T = "T cells CD8",
  CD4T = c(
    "T cells CD4 naive",
    "T cells CD4 memory resting",
    "T cells CD4 memory activated",
    "T cells follicular helper",
    "T cells regulatory (Tregs)"
  ),
  B = c("B cells naive", "B cells memory", "Plasma cells"),
  NK = c("NK cells resting", "NK cells activated")
)

state_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_state_map[[x]])
  names(out) <- target
  out
}
citmic_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_citmic_map[[x]])
  names(out) <- target
  out
}
immucellai_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_immucellai_map[[x]])
  names(out) <- target
  out
}
lm22_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_lm22_map[[x]])
  names(out) <- target
  out
}

read_external_matrix <- function(file) {
  first <- readLines(file, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  dat <- read.table(
    file, sep = sep, header = TRUE, row.names = 1L, check.names = FALSE,
    quote = "", comment.char = ""
  )
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  mat
}

evaluate_prediction <- function(pred.mat, truth.mat, target.mapper) {
  pred.mat <- as.matrix(pred.mat)
  if (!all(rownames(truth.mat) %in% rownames(pred.mat)) &&
      all(rownames(truth.mat) %in% colnames(pred.mat))) {
    pred.mat <- t(pred.mat)
  }
  pred.mat <- pred.mat[rownames(truth.mat), , drop = FALSE]
  rows <- lapply(colnames(truth.mat), function(ct) {
    states <- intersect(target.mapper(ct)[[ct]], colnames(pred.mat))
    truth <- truth.mat[, ct]
    pred <- if (length(states) == 0) rep(NA_real_, nrow(truth.mat)) else
      rowSums(pred.mat[, states, drop = FALSE], na.rm = TRUE)
    ok <- is.finite(pred) & is.finite(truth)
    fit <- if (sum(ok) >= 3 && stats::sd(pred[ok]) > 0 && stats::sd(truth[ok]) > 0) {
      stats::coef(stats::lm(pred[ok] ~ truth[ok]))
    } else {
      c(`(Intercept)` = NA_real_, truth = NA_real_)
    }
    data.frame(
      CellType = ct,
      PredictedStates = paste(states, collapse = "+"),
      Pearson = safe_cor(pred, truth, "pearson"),
      Spearman = safe_cor(pred, truth, "spearman"),
      RMSE = if (sum(ok) > 0) sqrt(mean((pred[ok] - truth[ok])^2)) else NA_real_,
      MAE = if (sum(ok) > 0) mean(abs(pred[ok] - truth[ok])) else NA_real_,
      Slope = unname(fit[2]),
      Intercept = unname(fit[1]),
      N = sum(ok),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

message("Reading matrices...")
bulk_all <- read_external_matrix(bulk_file)
reference_all <- read_external_matrix(reference_file)
truth_all <- read_external_matrix(truth_file)
bulk_all[!is.finite(bulk_all) | bulk_all < 0] <- 0
reference_all[!is.finite(reference_all) | reference_all < 0] <- 0
if (!is.na(sample_limit)) {
  bulk_all <- bulk_all[, seq_len(min(sample_limit, ncol(bulk_all))), drop = FALSE]
}

common <- intersect(rownames(bulk_all), rownames(reference_all))
bulk_all <- bulk_all[common, , drop = FALSE]
reference_all <- reference_all[common, , drop = FALSE]
keep <- rowSums(bulk_all, na.rm = TRUE) > 0 & rowSums(reference_all, na.rm = TRUE) > 0
bulk_all <- bulk_all[keep, , drop = FALSE]
reference_all <- reference_all[keep, , drop = FALSE]

common_samples <- intersect(colnames(bulk_all), rownames(truth_all))
if (length(common_samples) < 3L) stop("Fewer than three aligned bulk/truth samples.")
bulk_all <- bulk_all[, common_samples, drop = FALSE]
truth_all <- truth_all[common_samples, , drop = FALSE]
truth_all <- truth_all[, intersect(colnames(truth_all), names(external_state_map)), drop = FALSE]
if (ncol(truth_all) == 0L) stop("No evaluable truth cell types.")

hierarchy_base <- create_default_53_hierarchy(colnames(reference_all), hierarchy.mode = "tcell")
hierarchy_base <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_base, reference_all)
reference_all <- reference_all[, hierarchy_base$state_name, drop = FALSE]

bulk_external <- bulk_all
reference_external <- reference_all
external_gene_strategy <- "all_common_genes"

immune_genes <- intersect(read_gene_list(immune_gene_file), rownames(reference_all))
if (length(immune_genes) >= 100) {
  bulk_immucellai2 <- bulk_all[immune_genes, , drop = FALSE]
  reference_immucellai2 <- reference_all[immune_genes, , drop = FALSE]
  immucellai2_gene_strategy <- "user_immune_genes"
} else {
  bulk_immucellai2 <- bulk_all
  reference_immucellai2 <- reference_all
  immucellai2_gene_strategy <- "all_common_genes"
}

message("Dataset: ", dataset_name,
        "; samples: ", ncol(bulk_all), "; targets: ", ncol(truth_all),
        "; external genes: ", nrow(bulk_external),
        "; ImmuCellAI2 genes: ", nrow(bulk_immucellai2),
        "; cores: ", n_cores,
        "; ImmuCellAI2 gene_strategy: ", immucellai2_gene_strategy,
        "; external gene_strategy: ", external_gene_strategy)

metric_list <- list()
summary_list <- list()
run_status <- data.frame(Method = character(), Status = character(), Message = character())

add_status <- function(method, status, msg = "") {
  run_status <<- rbind(run_status, data.frame(Method = method, Status = status, Message = msg))
  write.table(run_status, file.path(out_dir, "run_status.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

add_method <- function(method, pred.mat, mapper = state_target_mapper, genes.used = nrow(bulk_external)) {
  metrics <- evaluate_prediction(pred.mat, truth_all, mapper)
  metrics$Method <- method
  metric_list[[method]] <<- metrics[, c("Method", setdiff(colnames(metrics), "Method"))]
  summary_list[[method]] <<- summarize_eval(metrics, method, genes.used)
  write.table(pred.mat, file.path(out_dir, paste0(method, "_fraction.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(metric_list[[method]], file.path(out_dir, paste0(method, "_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  add_status(method, "OK", "")
}

safe_run <- function(method, expr) {
  dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)
  dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
  add_status(method, "RUNNING", "")
  tryCatch(expr, error = function(e) {
    add_status(method, "FAILED", conditionMessage(e))
    message("FAILED ", method, ": ", conditionMessage(e))
  })
}

reuse_method <- function(method, mapper = state_target_mapper, genes.used = nrow(bulk_external)) {
  fraction_file <- file.path(out_dir, paste0(method, "_fraction.txt"))
  if (!file.exists(fraction_file)) return(FALSE)
  message("Reusing existing result for ", method)
  pred <- read.delim(fraction_file, row.names = 1, check.names = FALSE)
  add_method(method, pred, mapper = mapper, genes.used = genes.used)
  add_status(method, "REUSED", fraction_file)
  TRUE
}

if (isTRUE(run_immucellai2) &&
    !reuse_method(immucellai2_method_name, genes.used = nrow(bulk_immucellai2))) {
  safe_run(immucellai2_method_name, {
  message("Running ", immucellai2_method_name, "...")
  hierarchy <- create_default_53_hierarchy(colnames(reference_immucellai2), hierarchy.mode = immucellai2_mode)
  hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference_immucellai2)
  fit <- deconvolve_bulk_matrix(
    bulk.mat = bulk_immucellai2,
    reference = reference_immucellai2,
    hierarchy = hierarchy,
    hierarchy.mode = immucellai2_mode,
    inference.method = if (immucellai2_mode == "flat") "flat_vb" else "vb",
    add.unknown = immucellai2_add_unknown,
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    n.cores = n_cores,
    seed = 123
  )
  write_deconvolution_outputs(fit, file.path(out_dir, paste0(immucellai2_method_name, "_outputs")))
  add_method(immucellai2_method_name, fit$state.fraction, genes.used = nrow(bulk_immucellai2))
  })
}

if (isTRUE(run_bayesprism) &&
    !reuse_method("BayesPrism_first_state_chain600_burn500")) {
  safe_run("BayesPrism_first_state_chain600_burn500", {
  message("Running BayesPrism chain.length=", chain_length, ", burn.in=", burn_in, "...")
  hierarchy_bp <- create_default_53_hierarchy(colnames(reference_external), hierarchy.mode = "tcell")
  hierarchy_bp <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_bp, reference_external)
  reference.bp <- t(reference_external[, hierarchy_bp$state_name, drop = FALSE])
  mixture.bp <- t(bulk_external)
  set.seed(123)
  prism <- new.prism(
    reference = reference.bp,
    input.type = "count.matrix",
    cell.type.labels = hierarchy_bp$major_lineage,
    cell.state.labels = hierarchy_bp$state_name,
    key = NA_character_,
    mixture = mixture.bp,
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
  bp.state <- get.fraction(bp, which.theta = "first", state.or.type = "state")
  add_method("BayesPrism_first_state_chain600_burn500", bp.state)
  })
}

if (isTRUE(run_dwls) &&
    !reuse_method("DWLS_native_ref_top1000")) {
  safe_run("DWLS_native_ref_top1000", {
  message("Running native DWLS::solveDampenedWLS on top ", cibersort_signature_genes,
          " reference genes...")
  add_method("DWLS_native_ref_top1000", run_dwls_matrix(bulk_external, reference_external))
  })
}

if (isTRUE(run_music) &&
    !reuse_method("MuSiC_native_music.basic")) {
  safe_run("MuSiC_native_music.basic", {
  message("Running native MuSiC::music.basic on the 53-cell reference matrix...")
  add_method("MuSiC_native_music.basic", run_music_basic_matrix(bulk_external, reference_external))
  })
}

if (isTRUE(run_citmic) &&
    !reuse_method("CITMIC_native", citmic_target_mapper)) {
  safe_run("CITMIC_native", {
  message("Running CITMIC native with built-in signatures...")
  add_status("CITMIC_native_reference_note", "INFO", "CITMIC native API does not accept the 53-cell reference matrix; using package default signatures.")
  add_method("CITMIC_native", run_citmic_native_matrix(bulk_external), citmic_target_mapper)
  })
}

if (isTRUE(run_cibersort) &&
    !reuse_method("CIBERSORT_default_LM22", lm22_target_mapper, genes.used = 547)) {
  safe_run("CIBERSORT_default_LM22", {
  message("Running native CIBERSORT::cibersort with built-in LM22 and default parameters...")
  add_method("CIBERSORT_default_LM22", run_cibersort_lm22_matrix(bulk_external),
             lm22_target_mapper, genes.used = 547)
  })
}

if (isTRUE(run_immucellai) &&
    !reuse_method("ImmuCellAI_native", immucellai_target_mapper)) {
  safe_run("ImmuCellAI_native", {
  message("Running ImmuCellAI native with built-in signatures...")
  add_status("ImmuCellAI_native_reference_note", "INFO", "ImmuCellAI_new native API does not accept the 53-cell reference matrix; using package default signatures.")
  add_method("ImmuCellAI_native", run_immucellai_native_matrix(bulk_external), immucellai_target_mapper)
  })
}

if (length(metric_list) == 0) stop("No method completed successfully.")

all_metrics <- do.call(rbind, metric_list)
write.table(all_metrics, file.path(out_dir, "all_methods_per_celltype_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- do.call(rbind, summary_list)
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

print(summary_table)
message("Done. Results written to: ", out_dir)
message("Summary: ", file.path(out_dir, "comparison_summary_7tools.txt"))
message("Per-celltype Pearson wide table: ", file.path(out_dir, "per_celltype_pearson_7tools_wide.txt"))
