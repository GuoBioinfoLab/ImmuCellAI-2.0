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
  get <- function(path) read_tsv(file.path(data_dir, path))
  index <- get("panel_index.tsv")
  scores <- get("shared/Fig4_H_ROC_scores.tsv")
  predictions <- get("shared/Fig4_H_predictions.tsv")
  stopifnot(identical(index$panel, LETTERS[1:8]), nrow(scores) == 334L,
            !anyDuplicated(scores$sample), sum(scores$ResponseBinary) == 120L,
            setequal(predictions$sample, scores$sample), setequal(predictions$fold, 1:5))
  check_num(scores$Fine53_RF, predictions$Unified_ICB_Score[match(scores$sample, predictions$sample)], "Frozen OOF scores")
  list(index = index, scores = scores, predictions = predictions,
       verification = list(n = 334L, responders = 120L, panels = 8L, model_refitted = FALSE))
}
analyze_panel <- function(panel) {
  get <- function(path) read_tsv(file.path(data_dir, path))
  meta <- context$index[context$index$panel == panel, , drop = FALSE]
  scores <- get(paste0(panel, "/ROC_scores.tsv"))
  curve <- get(paste0(panel, "/ROC_plot_data.tsv"))
  expected <- context$scores[context$scores$SRA_study == meta$study, , drop = FALSE]
  check_table(scores, expected, "sample", c("SRA_study", "ResponseBinary", "Fine53_RF"))
  check_num(scores$fold, context$predictions$fold[match(scores$sample, context$predictions$sample)], "OOF folds")
  stopifnot(nrow(scores) == meta$N, sum(scores$ResponseBinary) == meta$responders)
  value <- auc_rank(scores$ResponseBinary, scores$Fine53_RF)
  check_num(value, meta$recomputed_AUC, "AUC")
  stopifnot(sprintf("%.3f", value) == sprintf("%.3f", meta$figure_AUC))
  ix <- order(-scores$Fine53_RF, method = "radix")
  y <- scores$ResponseBinary[ix]; s <- scores$Fine53_RF[ix]
  ends <- c(which(diff(s) != 0), length(s))
  expected_curve <- data.frame(threshold = c(Inf, s[ends], -Inf),
                               FPR = c(0, cumsum(y == 0)[ends] / sum(y == 0), 1),
                               TPR = c(0, cumsum(y == 1)[ends] / sum(y == 1), 1))
  for (column in names(expected_curve)) check_num(curve[[column]], expected_curve[[column]], paste(panel, column))
  vertices <- get(paste0(panel, "/source_pdf_vertices.tsv"))
  xx <- (vertices$x_pt - min(vertices$x_pt)) / diff(range(vertices$x_pt))
  yy <- (max(vertices$top_pt) - vertices$top_pt) / diff(range(vertices$top_pt))
  check_num(xx, head(curve$FPR, -1), "Historical PDF x", 1e-4)
  check_num(yy, head(curve$TPR, -1), "Historical PDF y", 1e-4)
  list(meta = meta, scores = scores, curve = curve, AUC = value)
}
