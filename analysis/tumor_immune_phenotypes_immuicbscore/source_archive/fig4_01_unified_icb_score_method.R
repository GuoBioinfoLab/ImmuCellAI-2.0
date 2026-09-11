suppressPackageStartupMessages({
  library(data.table)
})

# Unified ICB Score
# -----------------
# Definition used in the Fig4 analyses:
#   1. Start from an ImmuCellAI2 cell-fraction matrix with samples in rows and
#      the 53 immune cell types in columns.
#   2. Apply CLR transform: log(pmax(fraction, 0) + 1e-5) minus the row mean.
#   3. Standardize features using training-set mean and SD.
#   4. Fit a probability random forest with balanced class weights.
#   5. The predicted probability of response is the Unified ICB Score.

clr_transform <- function(x, eps = 1e-5) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx, na.rm = TRUE), "-")
}

class_weights <- function(y) {
  y <- as.integer(y)
  if (!all(y %in% c(0L, 1L))) stop("y must be encoded as 0/1.")
  if (length(unique(y)) < 2) stop("Both response classes are required.")
  ifelse(y == 1L, 0.5 / mean(y == 1L), 0.5 / mean(y == 0L))
}

fit_feature_scaler <- function(x_train) {
  center <- colMeans(x_train, na.rm = TRUE)
  scalev <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(center = center, scale = scalev)
}

apply_feature_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

train_unified_icb_score <- function(cell_fraction,
                                    response,
                                    n_trees = 1200,
                                    min_node_size = 8,
                                    n_threads = 8,
                                    seed = 123) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required to train Unified ICB Score.")
  }
  y <- as.integer(response)
  if (!all(y %in% c(0L, 1L))) stop("response must be encoded as 0/1.")
  x_clr <- clr_transform(cell_fraction, eps = 1e-5)
  scaler <- fit_feature_scaler(x_clr)
  x_scaled <- apply_feature_scaler(x_clr, scaler)
  cw <- class_weights(y)

  set.seed(seed)
  fit <- ranger::ranger(
    x = as.data.frame(x_scaled, check.names = FALSE),
    y = factor(y, levels = c(0, 1)),
    probability = TRUE,
    classification = TRUE,
    num.trees = n_trees,
    mtry = max(1, floor(sqrt(ncol(x_scaled)))),
    min.node.size = min_node_size,
    class.weights = c("0" = cw[which(y == 0L)[1]], "1" = cw[which(y == 1L)[1]]),
    seed = seed,
    num.threads = n_threads
  )

  structure(
    list(
      model = fit,
      celltypes = colnames(x_clr),
      scaler = scaler,
      eps = 1e-5,
      n_trees = n_trees,
      min_node_size = min_node_size
    ),
    class = "UnifiedICBScoreModel"
  )
}

predict_unified_icb_score <- function(object, cell_fraction) {
  if (!inherits(object, "UnifiedICBScoreModel")) {
    stop("object must be returned by train_unified_icb_score().")
  }
  missing_cells <- setdiff(object$celltypes, colnames(cell_fraction))
  if (length(missing_cells) > 0) {
    stop("Missing required cell types: ", paste(missing_cells, collapse = ", "))
  }
  x <- cell_fraction[, object$celltypes, drop = FALSE]
  x_clr <- clr_transform(x, eps = object$eps)
  x_scaled <- apply_feature_scaler(x_clr, object$scaler)
  as.numeric(predict(object$model, data = as.data.frame(x_scaled, check.names = FALSE))$predictions[, "1"])
}

write_unified_icb_score_method_note <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  note <- c(
    "Unified ICB Score method",
    "",
    "Input: ImmuCellAI2 53-cell fraction matrix, samples in rows and cell types in columns.",
    "Feature transform: CLR = log(pmax(fraction, 0) + 1e-5) - row mean of log-transformed fractions.",
    "Model: probability random forest implemented by ranger.",
    "Class balance: class.weights = 0.5 / class prevalence for responders and non-responders.",
    "Tuning used for the selected model: num.trees = 1200, mtry = floor(sqrt(53)), min.node.size = 8.",
    "Output score: predicted probability of response, named Unified ICB Score.",
    "",
    "In Fig4, model selection identified the selected broad model as:",
    "  feature_set = cell_clr53",
    "  method = ranger_balanced",
    "  model id = cell_clr53__ranger_balanced"
  )
  writeLines(note, file.path(out_dir, "Unified_ICB_score_method_note.txt"), useBytes = TRUE)
}

if (identical(environment(), globalenv()) && !interactive()) {
  out_dir <- "<LOCAL_R_ROOT>/Fig4/Unified_ICB_score_cluster_ICB_integration"
  write_unified_icb_score_method_note(out_dir)
  message("Unified ICB Score method note written to: ", out_dir)
}
