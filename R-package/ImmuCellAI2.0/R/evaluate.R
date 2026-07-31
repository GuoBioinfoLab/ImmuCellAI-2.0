#' Evaluate increasing-proportion simulation
#'
#' @param state.mat Sample x state prediction matrix.
#' @param truth.mat Optional state x sample truth matrix.
#' @param sample.info data.frame with Sample, TargetCell, and either Truth or TargetProportion.
#' @param intersection.only If TRUE, only evaluate targets directly present in prediction columns.
#' @return list with sample.info, metrics, matched.targets, unmatched.targets.
#' @export
evaluate_increasing_simulation <- function(state.mat, truth.mat = NULL, sample.info,
                                           intersection.only = TRUE) {
  state.mat <- as.matrix(state.mat)
  if (!all(c("Sample", "TargetCell") %in% colnames(sample.info))) {
    stop("sample.info must contain Sample and TargetCell columns.")
  }
  if (!"Truth" %in% colnames(sample.info)) {
    if ("TargetProportion" %in% colnames(sample.info)) {
      sample.info$Truth <- sample.info$TargetProportion
    } else if (!is.null(truth.mat)) {
      truth.mat <- as.matrix(truth.mat)
      sample.info$Truth <- mapply(function(ct, smp) {
        if (ct %in% rownames(truth.mat) && smp %in% colnames(truth.mat)) truth.mat[ct, smp] else NA_real_
      }, sample.info$TargetCell, sample.info$Sample)
    } else {
      stop("Provide sample.info$Truth, sample.info$TargetProportion, or truth.mat.")
    }
  }

  simulated.target.cells <- unique(sample.info$TargetCell)
  matched.targets <- intersect(simulated.target.cells, colnames(state.mat))
  unmatched.targets <- setdiff(simulated.target.cells, colnames(state.mat))

  sample.info$Predicted <- NA_real_
  sample.info$IsMatchedTarget <- sample.info$TargetCell %in% matched.targets
  for (i in seq_len(nrow(sample.info))) {
    smp <- sample.info$Sample[i]
    target <- sample.info$TargetCell[i]
    if (!smp %in% rownames(state.mat)) next
    if (!target %in% colnames(state.mat)) next
    sample.info$Predicted[i] <- state.mat[smp, target]
  }

  calc_metrics <- function(df) {
    truth <- df$Truth
    pred <- df$Predicted
    keep <- !is.na(truth) & !is.na(pred)
    truth <- truth[keep]
    pred <- pred[keep]
    if (length(truth) < 3) {
      return(data.frame(N = length(truth), Correlation = NA_real_, RMSE = NA_real_,
                        MAE = NA_real_, Bias = NA_real_, Slope = NA_real_, Intercept = NA_real_))
    }
    corr <- if (sd(truth) == 0 || sd(pred) == 0) NA_real_ else stats::cor(truth, pred)
    rmse <- sqrt(mean((truth - pred)^2))
    mae <- mean(abs(truth - pred))
    bias <- mean(pred - truth)
    lm.fit <- stats::lm(pred ~ truth)
    data.frame(N = length(truth), Correlation = corr, RMSE = rmse, MAE = mae,
               Bias = bias, Slope = stats::coef(lm.fit)[2], Intercept = stats::coef(lm.fit)[1])
  }

  targets.to.eval <- if (intersection.only) matched.targets else simulated.target.cells
  metrics.list <- list()
  for (ct in targets.to.eval) {
    df <- sample.info[sample.info$TargetCell == ct, , drop = FALSE]
    m <- calc_metrics(df)
    m$TargetCell <- ct
    metrics.list[[ct]] <- m
  }
  metrics <- if (length(metrics.list) > 0) do.call(rbind, metrics.list) else data.frame()
  if (nrow(metrics) > 0) {
    metrics <- metrics[, c("TargetCell", "N", "Correlation", "RMSE", "MAE", "Bias", "Slope", "Intercept")]
    metrics <- metrics[order(metrics$Correlation, decreasing = TRUE), ]
  }

  list(sample.info = sample.info, metrics = metrics,
       matched.targets = matched.targets, unmatched.targets = unmatched.targets)
}
