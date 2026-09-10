# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

root <- file.path(fig3_work_dir(), "Fig2/Newman_Monaco_7tools_validation")
lm22_dir <- file.path(root, "CIBERSORT_default_LM22")

datasets <- list(
  RNA_Sieve_newman_pbmcs = list(
    lm22_fraction = "Newman_PBMC_CIBERSORT_default_LM22_fraction.txt",
    lm22_metrics = "Newman_PBMC_CIBERSORT_default_LM22_metrics.txt"
  ),
  RNA_Sieve_monaco_pbmcs = list(
    lm22_fraction = "Monaco_PBMC_CIBERSORT_default_LM22_fraction.txt",
    lm22_metrics = "Monaco_PBMC_CIBERSORT_default_LM22_metrics.txt"
  )
)

summarize_eval <- function(metrics, method, genes.used) {
  data.frame(
    Method = method,
    Level = "external_ground_truth",
    N_targets = nrow(metrics),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MedianSpearman = median(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MedianRMSE = median(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE),
    MedianMAE = median(metrics$MAE, na.rm = TRUE),
    MeanSlope = mean(metrics$Slope, na.rm = TRUE),
    MedianSlope = median(metrics$Slope, na.rm = TRUE),
    GenesUsed = genes.used,
    ValidPearsonTargets = sum(!is.na(metrics$Pearson)),
    PenalizedMeanPearson = mean(ifelse(is.na(metrics$Pearson), -1, metrics$Pearson))
  )
}

for (dataset in names(datasets)) {
  cfg <- datasets[[dataset]]
  result_dir <- file.path(root, dataset, "seven_tools_compare")
  all_file <- file.path(result_dir, "all_methods_per_celltype_metrics.txt")
  summary_file <- file.path(result_dir, "comparison_summary_7tools.txt")

  all_metrics <- read.delim(all_file, check.names = FALSE)
  all_metrics <- all_metrics[
    !all_metrics$Method %in% c("CIBERSORT_native_ref_top1000", "CIBERSORT_default_LM22"),
    ,
    drop = FALSE
  ]

  lm22_metrics <- read.delim(file.path(lm22_dir, cfg$lm22_metrics), check.names = FALSE)
  lm22_metrics$Dataset <- NULL
  lm22_metrics$Method <- "CIBERSORT_default_LM22"
  lm22_metrics$N <- lm22_metrics$N_samples
  lm22_metrics$N_samples <- NULL
  lm22_metrics <- lm22_metrics[, colnames(all_metrics), drop = FALSE]
  all_metrics <- rbind(all_metrics, lm22_metrics)

  old_summary <- read.delim(summary_file, check.names = FALSE)
  old_summary <- old_summary[
    !old_summary$Method %in% c("CIBERSORT_native_ref_top1000", "CIBERSORT_default_LM22"),
    ,
    drop = FALSE
  ]
  lm22_summary <- summarize_eval(lm22_metrics, "CIBERSORT_default_LM22", 547)
  summary_new <- rbind(old_summary, lm22_summary)
  summary_new <- summary_new[order(summary_new$PenalizedMeanPearson, decreasing = TRUE), , drop = FALSE]

  write.table(all_metrics, all_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_new, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(lm22_metrics, file.path(result_dir, "CIBERSORT_default_LM22_metrics.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  fraction <- read.delim(file.path(lm22_dir, cfg$lm22_fraction),
                         row.names = 1, check.names = FALSE)
  write.table(fraction, file.path(result_dir, "CIBERSORT_default_LM22_fraction.txt"),
              sep = "\t", quote = FALSE, col.names = NA)

  pearson_wide <- reshape(
    all_metrics[, c("CellType", "Method", "Pearson")],
    idvar = "CellType", timevar = "Method", direction = "wide"
  )
  colnames(pearson_wide) <- sub("^Pearson\\.", "", colnames(pearson_wide))
  write.table(pearson_wide, file.path(result_dir, "per_celltype_pearson_7tools_wide.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  spearman_wide <- reshape(
    all_metrics[, c("CellType", "Method", "Spearman")],
    idvar = "CellType", timevar = "Method", direction = "wide"
  )
  colnames(spearman_wide) <- sub("^Spearman\\.", "", colnames(spearman_wide))
  write.table(spearman_wide, file.path(result_dir, "per_celltype_spearman_7tools_wide.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  status_file <- file.path(result_dir, "run_status.txt")
  if (file.exists(status_file)) {
    status <- read.delim(status_file, check.names = FALSE)
    status <- status[
      !status$Method %in% c("CIBERSORT_native_ref_top1000", "CIBERSORT_default_LM22"),
      ,
      drop = FALSE
    ]
    status <- rbind(
      status,
      data.frame(
        Method = "CIBERSORT_default_LM22",
        Status = "REPLACED_WITH_DEFAULT_LM22",
        Message = "Built-in LM22; default maxSize=500, perm=0, QN=TRUE"
      )
    )
    write.table(status, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

dataset_names <- names(datasets)
combined <- do.call(rbind, lapply(dataset_names, function(dataset) {
  x <- read.delim(file.path(root, dataset, "seven_tools_compare", "comparison_summary_7tools.txt"),
                  check.names = FALSE)
  x$Dataset <- dataset
  x
}))
combined <- combined[, c("Dataset", setdiff(colnames(combined), "Dataset"))]
write.table(combined, file.path(root, "Newman_Monaco_7tools_comparison_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

best <- do.call(rbind, lapply(split(combined, combined$Dataset), function(x) {
  x[which.max(x$PenalizedMeanPearson), , drop = FALSE]
}))
rownames(best) <- NULL
write.table(best, file.path(root, "Newman_Monaco_7tools_best_method.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("CIBERSORT results replaced with default LM22.\n\n")
print(combined[combined$Method %in% c(
  "ImmuCellAI2_flat_VB", "ImmuCellAI2_tcell_VB_UNKNOWN", "CIBERSORT_default_LM22"
), c("Dataset", "Method", "MeanPearson", "MedianPearson", "MeanRMSE", "MeanMAE")],
row.names = FALSE)
