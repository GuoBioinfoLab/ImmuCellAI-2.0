library(data.table)

cluster_dir <- "<LOCAL_R_ROOT>/Fig4/cluster_TCGA"
out_file <- file.path(cluster_dir, "cluster_TCGA_k_selection_metrics.txt")

read_class <- function(k) {
  f <- file.path(cluster_dir, sprintf("cluster_TCGA.k=%d.consensusClass.csv", k))
  x <- fread(f, header = FALSE, col.names = c("Sample", "Cluster"))
  x[, Cluster := as.integer(Cluster)]
  x
}

evaluate_k <- function(k) {
  message("Evaluating k=", k)
  class_dt <- read_class(k)
  mat_file <- file.path(cluster_dir, sprintf("cluster_TCGA.k=%d.consensusMatrix.csv", k))
  mat_dt <- fread(mat_file, check.names = FALSE)

  row_ids <- mat_dt[[1L]]
  mat_dt[[1L]] <- NULL
  mat <- as.matrix(mat_dt)
  storage.mode(mat) <- "double"
  rm(mat_dt)
  gc()

  n <- nrow(mat)
  if (length(row_ids) != nrow(class_dt)) {
    warning("k=", k, ": consensus matrix row count and class count differ.")
  }

  # ConsensusClusterPlus writes matrix rows in the same item order as class output
  # in this CSV export, but keep an explicit sanity check when row ids are sample IDs.
  if (all(row_ids %in% class_dt$Sample)) {
    class_ordered <- class_dt[match(row_ids, Sample)]
  } else {
    class_ordered <- class_dt
  }

  lower_idx <- lower.tri(mat)
  vals <- mat[lower_idx]
  cls <- class_ordered$Cluster
  same_cluster <- outer(cls, cls, "==")[lower_idx]

  sizes <- as.integer(table(factor(cls, levels = sort(unique(cls)))))
  names(sizes) <- sort(unique(cls))

  within_vals <- vals[same_cluster]
  between_vals <- vals[!same_cluster]

  data.frame(
    K = k,
    N = n,
    ClusterSizes = paste(sprintf("%s:%s", names(sizes), sizes), collapse = "; "),
    MinClusterSize = min(sizes),
    MaxClusterSize = max(sizes),
    SizeImbalance = max(sizes) / min(sizes),
    PAC_0.1_0.9 = mean(vals > 0.1 & vals < 0.9),
    PAC_0.05_0.95 = mean(vals > 0.05 & vals < 0.95),
    ConfidentLe0.1Ge0.9 = mean(vals <= 0.1 | vals >= 0.9),
    MeanConsensus = mean(vals),
    MedianConsensus = median(vals),
    WithinMean = mean(within_vals),
    WithinMedian = median(within_vals),
    BetweenMean = mean(between_vals),
    BetweenMedian = median(between_vals),
    Separation = mean(within_vals) - mean(between_vals),
    stringsAsFactors = FALSE
  )
}

ks <- 2:6
metrics <- rbindlist(lapply(ks, evaluate_k), fill = TRUE)
metrics[, DeltaPAC := c(NA_real_, diff(PAC_0.1_0.9))]
metrics[, DeltaSeparation := c(NA_real_, diff(Separation))]
fwrite(metrics, out_file, sep = "\t", quote = FALSE)
print(metrics)
cat("\nSaved:", out_file, "\n")
