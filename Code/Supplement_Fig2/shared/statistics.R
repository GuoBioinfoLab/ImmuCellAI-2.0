read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                     fileEncoding = "UTF-8-BOM", quote = "\"", comment.char = "", na.strings = c("", "NA", "NaN", "nan"))
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
    else if (!identical(as.character(a[[column]]), as.character(b[[column]]))) stop("Table column differs: ", column)
  }
  invisible(TRUE)
}
auc_rank <- function(y, score) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  stopifnot(n1 > 0, n0 > 0, all(is.finite(score)), all(y %in% c(0, 1)))
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

methods <- c("ImmuCellAI2", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC_basic")
datasets <- c("GSE164522", "GSE146771", "GSE176078", "GSE162025")
prepare_analysis <- function(data_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
  get <- function(path) read_tsv(file.path(data_dir, path))
  shown <- get("_shared/pie_plot_data.tsv"); full <- get("_shared/all_computed_correlations.tsv")
  keys <- c("dataset", "method", "target")
  stopifnot(nrow(shown) == 469L, nrow(full) == 490L, !anyDuplicated(table_key(shown, keys)), !anyDuplicated(table_key(full, keys)))
  expected <- full[!(full$dataset == "GSE164522" & full$target %in% c("Th1", "Th17", "NKT")), ]
  check_table(shown, expected, keys)
  old <- get("_shared/previous_four_cohort_plot_data.tsv")
  check_table(shown[shown$dataset != "GSE164522", ], old[old$dataset != "GSE164522", ], keys)
  list(shown = shown, full = full, verification = list(computed = 490L, displayed = 469L, inference_rerun = FALSE))
}
analyze_cohort <- function(panel) {
  dataset <- datasets[match(panel, LETTERS[1:4])]
  root <- file.path(data_dir, "_shared/cohorts", dataset); folder <- file.path(root, "seven_tools/full")
  protocol <- jsonlite::fromJSON(file.path(folder, "input/locked_protocol.json"))
  mapping_path <- file.path(folder, "input/locked_mapping.tsv")
  stopifnot(tolower(unname(tools::md5sum(mapping_path))) == tolower(protocol$mapping_md5))
  for (name in c("sample_manifest.tsv", "truth_evaluation_fractions.tsv", "simulation_protocol.json"))
    stopifnot(tolower(unname(tools::md5sum(file.path(root, name)))) == tolower(protocol$source_md5[[name]]))
  manifest <- read_tsv(file.path(root, "sample_manifest.tsv")); truth <- read_tsv(file.path(root, "truth_evaluation_fractions.tsv"))
  mapping <- read_tsv(mapping_path); metrics <- read_tsv(file.path(folder, "metrics/per_celltype_correlations.tsv"))
  obs <- read_tsv(file.path(folder, "metrics/sample_level_truth_prediction.tsv"))
  steps <- if (panel == "A") 40L else 41L
  targets <- c(29L, 21L, 9L, 11L)[match(panel, LETTERS[1:4])]
  stopifnot(nrow(manifest) == steps * targets, nrow(manifest) == protocol$n_samples, !anyDuplicated(manifest$sample_id),
            identical(manifest$sample_id, truth$sample_id), identical(manifest$sample_id, protocol$sample_ids),
            nrow(metrics) == targets * 7L, nrow(obs) == 7L * nrow(manifest),
            !anyDuplicated(table_key(obs, c("method", "target", "sample_id"))))
  check_table(context$full[context$full$dataset == dataset, ], metrics, c("dataset", "method", "target"))
  shown <- read_tsv(file.path(data_dir, panel, "pie_plot_data.tsv"))
  check_table(shown, context$shown[context$shown$dataset == dataset, ], c("dataset", "method", "target"))
  predictions <- list()
  for (method in methods) {
    pred <- read_tsv(file.path(folder, "predictions", paste0(method, ".tsv")))
    stopifnot(identical(pred$sample_id, manifest$sample_id), all(is.finite(as.matrix(pred[-1]))))
    metadata <- jsonlite::fromJSON(file.path(folder, "predictions", paste0(method, "_metadata.json")))
    stopifnot(metadata$protocol_md5 == unname(tools::md5sum(file.path(folder, "input/locked_protocol.json"))))
    predictions[[method]] <- pred
  }
  if (panel == "A") {
    stopifnot(protocol$parameters$ImmuCellAI2[["n.iter"]] == 50L, identical(protocol$parameters$ImmuCellAI2[["add.unknown"]], FALSE))
    expected <- read_tsv(file.path(root, "expected_ImmuCellAI2_predictions.tsv"))
    stopifnot(nrow(expected) == 1160L, ncol(expected) == 54L, identical(names(expected), names(predictions$ImmuCellAI2)))
    check_table(predictions$ImmuCellAI2, expected, "sample_id")
  }
  errors <- numeric()
  for (i in seq_len(nrow(metrics))) {
    row <- metrics[i, ]; ix <- which(manifest$target_cell_type == row$target)
    gradient <- seq(if (steps == 40L) .205 else .20, .40, by = .005)
    stopifnot(length(ix) == steps); check_num(manifest$step[ix], seq_len(steps))
    check_num(manifest$target_fraction[ix], gradient); check_num(truth[[row$target]][ix], gradient)
    pairs <- obs[obs$method == row$method & obs$target == row$target, ]; pairs <- pairs[order(pairs$step), ]
    stopifnot(identical(pairs$sample_id, manifest$sample_id[ix]))
    check_num(pairs$truth, truth[[row$target]][ix])
    mapped <- mapping[mapping$method == row$method & mapping$target == row$target, ]; stopifnot(nrow(mapped) == 1L)
    states <- mapped$predicted_states
    if (is.na(states) || !nzchar(states)) {
      status <- "unavailable"; r <- NA_real_; stopifnot(all(is.na(pairs$prediction)))
    } else {
      states <- strsplit(states, " | ", fixed = TRUE)[[1]]
      values <- rowSums(predictions[[row$method]][ix, states, drop = FALSE])
      check_num(values, pairs$prediction)
      status <- if (sd(values) == 0) "constant_prediction" else "valid"
      r <- if (status == "valid") cor(values, truth[[row$target]][ix]) else NA_real_
    }
    stopifnot(row$status == status); check_num(row$Pearson, r, paste(dataset, row$method, row$target))
    if (is.finite(r)) errors <- c(errors, abs(row$Pearson - r))
  }
  list(data = shown, cells = read_tsv(file.path(data_dir, panel, "cell_order.tsv")), dataset = dataset,
       n_samples = nrow(manifest), maximum_error = max(errors))
}
