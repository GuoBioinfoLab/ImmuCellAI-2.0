args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript benchmarks/evaluate_predictions.R",
      "prediction_file truth_file output_dir [mapping_file]"
    )
  )
}

prediction_file <- args[1]
truth_file <- args[2]
output_dir <- args[3]
mapping_file <- if (length(args) >= 4L) args[4] else NA_character_
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(file) {
  utils::read.delim(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

pred <- read_tab(prediction_file)
if (ncol(pred) < 2L) stop("Prediction file must contain sample names and at least one cell type.")
names(pred)[1] <- "Sample"
pred$Sample <- as.character(pred$Sample)
if (anyDuplicated(pred$Sample)) stop("Prediction sample names are not unique.")

truth <- read_tab(truth_file)
if (ncol(truth) < 2L) stop("Truth file must contain sample names and at least one target.")

long_truth <- if (all(c("Sample", "TargetCell", "Truth") %in% names(truth))) {
  truth[, c("Sample", "TargetCell", "Truth")]
} else {
  names(truth)[1] <- "Sample"
  do.call(
    rbind,
    lapply(names(truth)[-1], function(target) {
      data.frame(
        Sample = truth$Sample,
        TargetCell = target,
        Truth = truth[[target]],
        stringsAsFactors = FALSE
      )
    })
  )
}
long_truth$Sample <- as.character(long_truth$Sample)
long_truth$TargetCell <- as.character(long_truth$TargetCell)
long_truth$Truth <- suppressWarnings(as.numeric(long_truth$Truth))

targets <- unique(long_truth$TargetCell)
mapping <- if (!is.na(mapping_file)) {
  x <- read_tab(mapping_file)
  required <- c("TargetCell", "PredictionColumns")
  if (!all(required %in% names(x))) {
    stop("Mapping file must contain TargetCell and PredictionColumns columns.")
  }
  x[, required]
} else {
  data.frame(
    TargetCell = targets,
    PredictionColumns = targets,
    stringsAsFactors = FALSE
  )
}
mapping$TargetCell <- as.character(mapping$TargetCell)
mapping$PredictionColumns <- as.character(mapping$PredictionColumns)

aggregate_prediction <- function(target) {
  spec <- mapping$PredictionColumns[match(target, mapping$TargetCell)]
  if (length(spec) == 0L || is.na(spec) || !nzchar(trimws(spec))) {
    return(rep(NA_real_, nrow(pred)))
  }
  columns <- trimws(strsplit(spec, "\\+", fixed = FALSE)[[1]])
  if (!all(columns %in% names(pred))) {
    return(rep(NA_real_, nrow(pred)))
  }
  values <- lapply(pred[, columns, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))
  mat <- as.matrix(as.data.frame(values, check.names = FALSE))
  if (ncol(mat) == 1L) as.numeric(mat[, 1]) else rowSums(mat)
}

prediction_long <- do.call(
  rbind,
  lapply(targets, function(target) {
    data.frame(
      Sample = pred$Sample,
      TargetCell = target,
      Predicted = aggregate_prediction(target),
      stringsAsFactors = FALSE
    )
  })
)

aligned <- merge(
  long_truth,
  prediction_long,
  by = c("Sample", "TargetCell"),
  all.x = TRUE,
  sort = FALSE
)

metric_one <- function(dat) {
  ok <- is.finite(dat$Truth) & is.finite(dat$Predicted)
  n_total <- nrow(dat)
  n_valid <- sum(ok)
  truth <- dat$Truth[ok]
  estimate <- dat$Predicted[ok]

  if (all(is.na(dat$Predicted))) {
    status <- "unavailable"
  } else if (n_valid < 3L) {
    status <- "insufficient_samples"
  } else if (stats::sd(estimate) == 0) {
    status <- "constant_prediction"
  } else {
    status <- "valid"
  }

  pearson <- if (status == "valid" && stats::sd(truth) > 0) {
    stats::cor(truth, estimate, method = "pearson")
  } else {
    NA_real_
  }
  spearman <- if (status == "valid" && stats::sd(truth) > 0) {
    suppressWarnings(stats::cor(truth, estimate, method = "spearman"))
  } else {
    NA_real_
  }

  fit <- if (status == "valid" && stats::sd(truth) > 0) {
    stats::lm(estimate ~ truth)
  } else {
    NULL
  }

  data.frame(
    Status = status,
    N_total = n_total,
    N_valid = n_valid,
    Pearson = pearson,
    Spearman = spearman,
    RMSE = if (n_valid > 0L) sqrt(mean((estimate - truth)^2)) else NA_real_,
    MAE = if (n_valid > 0L) mean(abs(estimate - truth)) else NA_real_,
    Bias = if (n_valid > 0L) mean(estimate - truth) else NA_real_,
    Slope = if (!is.null(fit)) unname(stats::coef(fit)[2]) else NA_real_,
    Intercept = if (!is.null(fit)) unname(stats::coef(fit)[1]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(
  rbind,
  lapply(split(aligned, aligned$TargetCell), function(dat) {
    cbind(TargetCell = dat$TargetCell[1], metric_one(dat), stringsAsFactors = FALSE)
  })
)
rownames(metrics) <- NULL

valid <- metrics$Status == "valid"
summary <- data.frame(
  N_targets = nrow(metrics),
  N_valid_Pearson = sum(is.finite(metrics$Pearson)),
  N_constant_prediction = sum(metrics$Status == "constant_prediction"),
  N_unavailable = sum(metrics$Status == "unavailable"),
  N_insufficient_samples = sum(metrics$Status == "insufficient_samples"),
  MeanPearson = if (any(valid)) mean(metrics$Pearson[valid], na.rm = TRUE) else NA_real_,
  MedianPearson = if (any(valid)) stats::median(metrics$Pearson[valid], na.rm = TRUE) else NA_real_,
  MeanSpearman = if (any(valid)) mean(metrics$Spearman[valid], na.rm = TRUE) else NA_real_,
  MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
  MeanMAE = mean(metrics$MAE, na.rm = TRUE),
  MeanBias = mean(metrics$Bias, na.rm = TRUE),
  stringsAsFactors = FALSE
)

utils::write.table(
  aligned,
  file.path(output_dir, "aligned_sample_target_values.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  metrics,
  file.path(output_dir, "per_celltype_metrics.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  summary,
  file.path(output_dir, "summary_metrics.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(summary)
