#' Default adaptive UNKNOWN parameters
#'
#' These are the tuned gated-fallback parameters used in the current package.
#'
#' @return Named list of adaptive parameters.
#' @export
default_adaptive_params <- function() {
  list(
    lambda.ref = 5000,
    learning.rate = 0.05,
    residual.quantile = 0.50,
    min.unknown.theta = 0.03,
    min.residual.fraction = 0.02,
    max.effective.residual.genes = 150,
    strong.residual.fraction = 0.05,
    strong.max.effective.genes = 120,
    fallback.learning.rate = 0.01,
    fallback.lambda.ref = 10000
  )
}

#' Extract simulation truth from sample names
#'
#' Assumes sample names such as CD4Tn_1, CD4Tn_2, ..., where the numeric suffix
#' indexes target proportions seq(0.005, 0.2, by = 0.005).
#'
#' @param sample.names Character vector of sample names.
#' @param truth.seq Numeric vector of target proportions.
#' @return data.frame with Sample, TargetCell, Index and Truth.
#' @export
extract_sample_truth_from_names <- function(
  sample.names,
  truth.seq = seq(0.005, 0.2, by = 0.005)
) {
  sample.info <- data.frame(
    Sample = sample.names,
    TargetCell = sub("_[0-9]+$", "", sample.names),
    Index = suppressWarnings(as.numeric(sub("^.*_([0-9]+)$", "\\1", sample.names))),
    stringsAsFactors = FALSE
  )
  sample.info$Truth <- truth.seq[sample.info$Index]
  sample.info
}

#' Build BayesPrism-panel output from state/subtype matrices
#'
#' This helper maps our output to a common BayesPrism-style panel: CD4Tn,
#' CD8Tn, FOB, ILC3, TAM, Treg, cDC1, cDC2, CMonocyte, NMonocyte, pDC,
#' and Thelper.
#'
#' @param state.mat Sample x state matrix.
#' @return Sample x panel matrix.
#' @export
build_bayesprism_panel <- function(state.mat) {
  state.mat <- as.matrix(state.mat)
  panel <- c(
    "CD4Tn", "CD8Tn", "FOB", "ILC3", "TAM", "Treg",
    "cDC1", "cDC2", "CMonocyte", "NMonocyte", "pDC", "Thelper"
  )
  out <- matrix(NA_real_, nrow = nrow(state.mat), ncol = length(panel))
  rownames(out) <- rownames(state.mat)
  colnames(out) <- panel

  get_state <- function(name) {
    if (name %in% colnames(state.mat)) state.mat[, name] else rep(NA_real_, nrow(state.mat))
  }

  out[, "CD4Tn"] <- get_state("CD4Tn")
  out[, "CD8Tn"] <- get_state("CD8Tn")
  out[, "FOB"] <- get_state("FOB")
  out[, "ILC3"] <- get_state("ILC3")
  out[, "TAM"] <- get_state("TAM")
  out[, "Treg"] <- get_state("Treg")
  out[, "cDC1"] <- get_state("cDC1")
  out[, "cDC2"] <- get_state("cDC2")
  out[, "CMonocyte"] <- if ("CMonocyte" %in% colnames(state.mat)) get_state("CMonocyte") else get_state("cMo")
  out[, "NMonocyte"] <- if ("NMonocyte" %in% colnames(state.mat)) get_state("NMonocyte") else get_state("ncMo")
  out[, "pDC"] <- get_state("pDC")

  helper.states <- intersect(c("Th1", "Th2", "Th17", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17"), colnames(state.mat))
  if (length(helper.states) > 0) {
    out[, "Thelper"] <- rowSums(state.mat[, helper.states, drop = FALSE])
  }

  out
}

#' Map target state names to BayesPrism-panel names
#'
#' @param x Character vector of target cell names.
#' @return Character vector of mapped target names.
#' @export
map_target_to_bayesprism_panel <- function(x) {
  y <- x
  y[x %in% c("cMo", "CMonocyte")] <- "CMonocyte"
  y[x %in% c("ncMo", "NMonocyte")] <- "NMonocyte"
  y[x %in% c("Th1", "Th2", "Th17", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17")] <- "Thelper"
  y
}

#' Map simulation target names to reference state names
#'
#' @param target Character vector of target cell names from sample names.
#' @return Named list when target has length > 1, otherwise a character vector.
#' @export
map_simulation_target_to_states <- function(target) {
  mapper <- list(
    CD4Tmemory = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
    CD8Tmemory = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
    CMonocyte = "cMo",
    NMonocyte = "ncMo",
    GC_B = "BGC",
    Thelper = c("Th1", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17", "Th17", "Th2"),
    cytotoxicNK = "cNK",
    exhausted_T = c("Tex", "exhausted_T"),
    macrophage = c("M0", "M1", "M2"),
    mast_cell = "Mast cell",
    plasma = c("PB", "PC"),
    regulatoryNK = "NKreg"
  )
  one <- function(x) {
    if (x %in% names(mapper)) mapper[[x]] else x
  }
  if (length(target) == 1) return(one(target))
  out <- lapply(target, one)
  names(out) <- target
  out
}

#' Evaluate prediction matrix against increasing-proportion truth with optional mapping
#'
#' @param pred.mat Sample x cell-type prediction matrix.
#' @param sample.info data.frame with Sample, TargetCell and Truth.
#' @param target.mapper Optional function mapping TargetCell to prediction columns.
#' @return list with sample.info and metrics.
#' @export
evaluate_prediction_matrix <- function(
  pred.mat,
  sample.info,
  target.mapper = identity
) {
  pred.mat <- as.matrix(pred.mat)
  if (!all(c("Sample", "TargetCell", "Truth") %in% colnames(sample.info))) {
    stop("sample.info must contain Sample, TargetCell and Truth columns.")
  }

  mapped <- target.mapper(sample.info$TargetCell)
  if (!is.list(mapped) && nrow(sample.info) == 1 && length(mapped) != 1) {
    mapped <- list(mapped)
  }
  if (is.list(mapped)) {
    sample.info$EvalTarget <- vapply(mapped, paste, collapse = "+", FUN.VALUE = character(1))
  } else {
    sample.info$EvalTarget <- as.character(mapped)
  }
  sample.info$Predicted <- NA_real_

  for (i in seq_len(nrow(sample.info))) {
    smp <- sample.info$Sample[i]
    target <- if (is.list(mapped)) mapped[[i]] else mapped[i]
    target <- intersect(target, colnames(pred.mat))
    if (smp %in% rownames(pred.mat) && length(target) > 0) {
      sample.info$Predicted[i] <- sum(pred.mat[smp, target])
    }
  }

  targets <- unique(sample.info$EvalTarget)
  metrics.list <- list()

  for (ct in targets) {
    df <- sample.info[
      sample.info$EvalTarget == ct &
        !is.na(sample.info$Truth) &
        !is.na(sample.info$Predicted),
      , drop = FALSE
    ]
    if (nrow(df) < 3) next

    truth <- df$Truth
    pred <- df$Predicted
    pearson <- if (stats::sd(truth) == 0 || stats::sd(pred) == 0) NA_real_ else stats::cor(truth, pred, method = "pearson")
    spearman <- if (stats::sd(truth) == 0 || stats::sd(pred) == 0) NA_real_ else suppressWarnings(stats::cor(truth, pred, method = "spearman"))
    fit <- stats::lm(pred ~ truth)

    metrics.list[[ct]] <- data.frame(
      CellType = ct,
      N = nrow(df),
      Pearson = pearson,
      Spearman = spearman,
      RMSE = sqrt(mean((pred - truth)^2)),
      MAE = mean(abs(pred - truth)),
      Bias = mean(pred - truth),
      Slope = unname(stats::coef(fit)[2]),
      Intercept = unname(stats::coef(fit)[1]),
      stringsAsFactors = FALSE
    )
  }

  metrics <- if (length(metrics.list) == 0) data.frame() else do.call(rbind, metrics.list)
  if (nrow(metrics) > 0) {
    rownames(metrics) <- NULL
    metrics <- metrics[order(metrics$Pearson, decreasing = TRUE), ]
  }

  list(sample.info = sample.info, metrics = metrics)
}

#' Summarize a metric table
#'
#' @param metrics Metrics table from evaluate_prediction_matrix().
#' @param method Method name.
#' @param level Evaluation level.
#' @return One-row summary data.frame.
#' @export
summarize_metrics <- function(metrics, method, level = "state") {
  if (is.null(metrics) || nrow(metrics) == 0) {
    return(data.frame(
      Method = method, Level = level, N_targets = 0,
      MeanPearson = NA_real_, MedianPearson = NA_real_, MeanSpearman = NA_real_,
      MeanRMSE = NA_real_, MedianRMSE = NA_real_, MeanMAE = NA_real_,
      MeanSlope = NA_real_, MedianSlope = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Method = method,
    Level = level,
    N_targets = nrow(metrics),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = stats::median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MedianRMSE = stats::median(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE),
    MeanSlope = mean(metrics$Slope, na.rm = TRUE),
    MedianSlope = stats::median(metrics$Slope, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Run Experiment 1 benchmark for current method variants
#'
#' This runs the current package under controlled settings:
#' flat no UNKNOWN, flat + UNKNOWN, hierarchical no UNKNOWN, hierarchical + UNKNOWN,
#' and optionally hybrid no UNKNOWN / hybrid + UNKNOWN.
#'
#' @param bulk.mat Gene x sample bulk matrix.
#' @param reference Gene x state reference matrix.
#' @param out.dir Output directory.
#' @param sample.info Optional sample truth table. If NULL, extracted from sample names.
#' @param include.hybrid Whether to include hybrid variants.
#' @param n.iter Number of iterations. Default 150.
#' @param burn.in Burn-in iterations. Default 50.
#' @param n.cores Number of cores.
#' @param seed Random seed.
#' @param adaptive.params Adaptive UNKNOWN parameter list. Defaults to tuned parameters.
#' @return list containing model results, state metrics, BayesPrism-panel metrics and summaries.
#' @export
run_experiment1_benchmark <- function(
  bulk.mat,
  reference,
  out.dir = "experiment1_current_method_benchmark",
  sample.info = NULL,
  include.hybrid = FALSE,
  n.iter = 150,
  burn.in = 50,
  n.cores = 1,
  seed = 123,
  adaptive.params = default_adaptive_params()
) {
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

  if (is.null(sample.info)) {
    sample.info <- extract_sample_truth_from_names(colnames(bulk.mat))
  }

  methods <- data.frame(
    Method = c("flat_no_UNKNOWN", "flat_UNKNOWN", "hierarchical_no_UNKNOWN", "hierarchical_UNKNOWN"),
    hierarchy.mode = c("flat", "flat", "hierarchical", "hierarchical"),
    add.unknown = c(FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  if (isTRUE(include.hybrid)) {
    methods <- rbind(methods, data.frame(
      Method = c("hybrid_no_UNKNOWN", "hybrid_UNKNOWN"),
      hierarchy.mode = c("hybrid", "hybrid"),
      add.unknown = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    ))
  }

  result.list <- list()
  state.metric.list <- list()
  panel.metric.list <- list()
  summary.list <- list()

  for (i in seq_len(nrow(methods))) {
    method <- methods$Method[i]
    mode <- methods$hierarchy.mode[i]
    add.unknown <- methods$add.unknown[i]

    message("Running ", method, " ...")
    method.dir <- file.path(out.dir, method)
    if (!dir.exists(method.dir)) dir.create(method.dir, recursive = TRUE)

    hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = mode)

    res <- deconvolve_bulk_matrix(
      bulk.mat = bulk.mat,
      reference = reference,
      hierarchy = hierarchy,
      hierarchy.mode = mode,
      add.unknown = add.unknown,
      unknown.state = "UNKNOWN",
      n.iter = n.iter,
      burn.in = burn.in,
      adaptive.params = adaptive.params,
      n.cores = n.cores,
      seed = seed
    )

    write_deconvolution_outputs(res, method.dir)

    state.eval <- evaluate_prediction_matrix(
      pred.mat = res$state.fraction,
      sample.info = sample.info,
      target.mapper = identity
    )

    panel.mat <- build_bayesprism_panel(res$state.fraction)
    panel.eval <- evaluate_prediction_matrix(
      pred.mat = panel.mat,
      sample.info = sample.info,
      target.mapper = map_target_to_bayesprism_panel
    )

    utils::write.table(state.eval$metrics, file.path(method.dir, "state_level_metrics.txt"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(panel.eval$metrics, file.path(method.dir, "bayesprism_panel_metrics.txt"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(panel.mat, file.path(method.dir, "bayesprism_panel_prediction.txt"),
                       sep = "\t", quote = FALSE, col.names = NA)

    result.list[[method]] <- res
    state.metric.list[[method]] <- state.eval$metrics
    panel.metric.list[[method]] <- panel.eval$metrics
    summary.list[[paste0(method, "_state")]] <- summarize_metrics(state.eval$metrics, method, "state")
    summary.list[[paste0(method, "_bayesprism_panel")]] <- summarize_metrics(panel.eval$metrics, method, "bayesprism_panel")
  }

  summary.table <- do.call(rbind, summary.list)
  rownames(summary.table) <- NULL
  utils::write.table(summary.table, file.path(out.dir, "experiment1_summary.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    results = result.list,
    state.metrics = state.metric.list,
    bayesprism.panel.metrics = panel.metric.list,
    summary = summary.table,
    sample.info = sample.info,
    methods = methods
  ))
}

#' Run a fast grid search for hierarchical marker/simplex methods
#'
#' This is the recommended first pass for the GSE146771 increasing-proportion
#' simulation because it can scan many hierarchical settings quickly before any
#' heavier Bayesian/MCMC comparison is run.
#'
#' @param bulk.mat Gene x sample bulk matrix.
#' @param reference Gene x state reference matrix.
#' @param out.dir Output directory.
#' @param marker.grid Marker counts to test.
#' @param include.simplex Include projected-gradient simplex methods.
#' @param include.marker.score Include marker-score methods.
#' @param n.iter Projected-gradient iterations for simplex methods.
#' @param sample.info Optional sample truth table.
#' @return list with summary and per-method metrics.
#' @export
run_fast_hierarchy_grid <- function(
  bulk.mat,
  reference,
  out.dir = "fast_hierarchy_grid",
  marker.grid = c(20, 30, 50, 100, 200, 400),
  include.simplex = TRUE,
  include.marker.score = TRUE,
  n.iter = 80,
  sample.info = NULL
) {
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  if (is.null(sample.info)) sample.info <- extract_sample_truth_from_names(colnames(bulk.mat))

  metric.list <- list()
  summary.list <- list()

  add_result <- function(method, pred.mat) {
    ev <- evaluate_prediction_matrix(
      pred.mat = pred.mat,
      sample.info = sample.info,
      target.mapper = map_simulation_target_to_states
    )
    metric.list[[method]] <<- ev$metrics
    summary <- summarize_metrics(ev$metrics, method, "simulation_target")
    summary$ValidPearsonTargets <- sum(!is.na(ev$metrics$Pearson))
    summary$PenalizedMeanPearson <- mean(ifelse(is.na(ev$metrics$Pearson), -1, ev$metrics$Pearson))
    summary.list[[method]] <<- summary
    utils::write.table(ev$metrics, file.path(out.dir, paste0(method, "_metrics.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  for (marker.n in marker.grid) {
    if (isTRUE(include.marker.score)) {
      message("Running marker_score_raw_", marker.n, " ...")
      res.score <- deconvolve_marker_score(
        bulk.mat = bulk.mat,
        reference = reference,
        hierarchy = create_default_53_hierarchy(colnames(reference), hierarchy.mode = "hybrid"),
        marker.n = marker.n,
        normalize.scores = FALSE
      )
      add_result(paste0("marker_score_raw_", marker.n), res.score$state.score)

      message("Running marker_score_fraction_", marker.n, " ...")
      res.frac <- deconvolve_marker_score(
        bulk.mat = bulk.mat,
        reference = reference,
        hierarchy = create_default_53_hierarchy(colnames(reference), hierarchy.mode = "hybrid"),
        marker.n = marker.n,
        normalize.scores = TRUE
      )
      add_result(paste0("marker_score_fraction_", marker.n), res.frac$state.fraction)
    }

    if (isTRUE(include.simplex)) {
      for (mode in c("flat", "hierarchical", "hybrid")) {
        message("Running simplex_", mode, "_", marker.n, " ...")
        res <- deconvolve_bulk_fast(
          bulk.mat = bulk.mat,
          reference = reference,
          hierarchy.mode = mode,
          marker.n = marker.n,
          n.iter = n.iter
        )
        add_result(paste0("simplex_", mode, "_", marker.n), res$state.fraction)
      }
    }
  }

  summary.table <- do.call(rbind, summary.list)
  rownames(summary.table) <- NULL
  summary.table <- summary.table[order(summary.table$PenalizedMeanPearson, decreasing = TRUE), ]
  utils::write.table(summary.table, file.path(out.dir, "fast_hierarchy_grid_summary.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(summary = summary.table, metrics = metric.list, sample.info = sample.info))
}
