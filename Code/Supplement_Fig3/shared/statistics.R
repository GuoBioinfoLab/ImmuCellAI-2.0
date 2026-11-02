read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                     fileEncoding = "UTF-8-BOM", quote = "\"", comment.char = "")
table_key <- function(x, cols) do.call(paste, c(x[cols], sep = "|"))
check_num <- function(a, b, label = "Numeric comparison", tolerance = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) stop(label, ": missingness or length differs")
  finite <- is.finite(a) & is.finite(b)
  if (any(is.infinite(a) != is.infinite(b)) || any(a[!finite & !is.na(a)] != b[!finite & !is.na(b)]))
    stop(label, ": nonfinite values differ")
  if (any(abs(a[finite] - b[finite]) > tolerance)) stop(label, ": max difference ", max(abs(a[finite] - b[finite])))
  invisible(TRUE)
}
check_table <- function(a, b, keys, cols = intersect(names(a), names(b)), tolerance = 1e-12) {
  ka <- table_key(a, keys); kb <- table_key(b, keys)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(kb), setequal(ka, kb))
  b <- b[match(ka, kb), , drop = FALSE]
  for (column in cols) {
    if (is.numeric(a[[column]]) && is.numeric(b[[column]])) check_num(a[[column]], b[[column]], column, tolerance)
    else stopifnot(identical(as.character(a[[column]]), as.character(b[[column]])))
  }
  invisible(TRUE)
}
auc_rank <- function(y, score) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  stopifnot(n1 > 0, n0 > 0, all(is.finite(score)), all(y %in% c(0, 1)))
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

prepare_analysis <- function(data_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install the R package jsonlite.")
  root <- file.path(data_dir, "_shared/source")
  get <- function(path) read_tsv(file.path(root, path))
  protocol <- jsonlite::fromJSON(file.path(root, "locked_protocol.json"))
  truth <- get("input/truth_published.tsv"); manifest <- get("input/sample_manifest.tsv")
  stopifnot(identical(truth[[1]], manifest[[1]]), nrow(truth) == 45L,
            !anyDuplicated(truth[[1]]), nrow(get("input/analysis_genes.tsv")) == 5418L)
  primary <- which(tolower(as.character(manifest$primary_complete_case)) == "true")
  stopifnot(length(primary) == 28L)
  groups <- protocol$groups
  actual <- sapply(groups, function(g) truth[[g$flow]])
  scenarios <- list(B = "B", DC = "DC", Monocyte = "Monocyte", TNK = c("T", "NK"))
  summary <- read_tsv(file.path(data_dir, "B/known_cell_MAE_summary.tsv"))
  unknown <- get("results/unknown_fraction_metrics.tsv")
  check_table(summary, get("results/known_cell_summary.tsv"), "job")
  sample_mae <- read_tsv(file.path(data_dir, "B/per_sample_MAE.tsv"))
  errors <- read_tsv(file.path(data_dir, "B/per_sample_group_errors.tsv"))
  scatter <- read_tsv(file.path(data_dir, "A/UNKNOWN_vs_omitted_fraction.tsv"))
  changes <- read_tsv(file.path(data_dir, "B/paired_MAE_changes.tsv"))
  check_table(changes, get("results/paired_MAE_changes.tsv"), "scenario")
  stopifnot(nrow(errors) == 1512L, nrow(sample_mae) == 224L, nrow(scatter) == 112L)
  for (scenario in names(scenarios)) {
    removed <- scenarios[[scenario]]; known <- setdiff(names(groups), removed)
    t <- rowSums(actual[primary, removed, drop = FALSE])
    maes <- list()
    for (on in c(FALSE, TRUE)) {
      job <- paste0(scenario, if (on) "_on" else "_off")
      pred <- get(paste0("results/", job, "/prediction.tsv"))
      stopifnot(identical(pred[[1]], truth[[1]]), all(is.finite(as.matrix(pred[-1]))),
                all(as.matrix(pred[-1]) >= 0), max(abs(rowSums(pred[-1]) - 1)) < 1e-9,
                ("UNKNOWN" %in% names(pred)) == on,
                !any(unlist(lapply(groups[removed], `[[`, "states")) %in% names(pred)))
      settings <- jsonlite::fromJSON(file.path(root, "results", job, "settings.json"))
      for (key in names(protocol$settings)) stopifnot(identical(settings[[key]], protocol$settings[[key]]))
      stopifnot(identical(settings[["add.unknown"]], on))
      estimate <- sapply(groups[known], function(g) rowSums(pred[g$states]))[primary, , drop = FALSE]
      err <- abs(estimate - actual[primary, known, drop = FALSE]) * 100
      mae <- rowMeans(err); maes[[if (on) "on" else "off"]] <- mae
      z <- sample_mae[sample_mae$job == job, ]
      stopifnot(identical(z$sample_id, truth[[1]][primary]), all(z$n_known_groups == length(known)))
      check_num(mae, z$mae_pp, paste(job, "sample MAE"))
      for (g in known) {
        z <- errors[errors$job == job & errors$group == g, ]
        stopifnot(identical(z$sample_id, truth[[1]][primary]))
        check_num(z$truth_fraction, actual[primary, g], "Group truth")
        check_num(z$prediction_fraction, estimate[, g], "Group prediction")
        check_num(z$absolute_error_pp, err[, g], "Group absolute error")
      }
      z <- summary[summary$job == job, ]
      rs <- vapply(known, function(g) cor(estimate[, g], actual[primary, g]), numeric(1))
      check_num(c(z$mae_pp, z$mean_pearson, z$sd_pearson), c(mean(mae), mean(rs), sd(rs)), "Known summary")
      y <- if (on) pred$UNKNOWN[primary] else rep(0, 28)
      slope <- cov(t, y) / var(t); intercept <- mean(y) - slope * mean(t)
      r <- if (on) cor(t, y) else NA_real_
      z <- unknown[unknown$job == job, ]
      check_num(unlist(z[c("r", "slope", "intercept", "mae", "mean_truth", "mean_pred")]),
                c(r, slope, intercept, mean(abs(t - y)), mean(t), mean(y)), "Unknown summary")
      if (on) {
        z <- scatter[scatter$scenario == scenario, ]
        stopifnot(identical(z$sample_id, truth[[1]][primary]))
        check_num(z$omitted_truth_fraction, t); check_num(z$UNKNOWN_weight_fraction, y)
        check_num(z$omitted_truth_percent, 100 * t); check_num(z$UNKNOWN_weight_percent, 100 * y)
      }
    }
    delta <- maes$off - maes$on
    z <- changes[changes$scenario == scenario, ]
    check_num(z$mae_reduction_pp, mean(delta)); stopifnot(z$improved_samples == sum(delta > 0))
  }
  metrics <- read_tsv(file.path(data_dir, "A/unknown_fraction_metrics.tsv"))
  check_table(metrics, unknown[unknown$job %in% paste0(names(scenarios), "_on"), ], "job")
  list(scatter = scatter, metrics = metrics, summary = summary,
       verification = list(samples = 28L, jobs = 8L, scatter_points = 112L, bars = 8L,
         bootstrap_CI = "Frozen original 4000-resample intervals checked against source; random draws not replayed."))
}
