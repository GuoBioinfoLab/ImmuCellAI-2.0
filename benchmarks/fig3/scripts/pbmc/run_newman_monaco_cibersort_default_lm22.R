# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(CIBERSORT))

root <- file.path(fig3_work_dir(), "Fig2/Newman_Monaco_7tools_validation")
out_dir <- file.path(root, "CIBERSORT_default_LM22")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  Newman_PBMC = list(
    expression = file.path(root, "RNA_Sieve_newman_pbmcs/standardized/expression.csv"),
    truth = file.path(root, "RNA_Sieve_newman_pbmcs/standardized/true_proportions.csv"),
    previous_metrics = file.path(
      root,
      "RNA_Sieve_newman_pbmcs/seven_tools_compare/CIBERSORT_native_ref_top1000_metrics.txt"
    )
  ),
  Monaco_PBMC = list(
    expression = file.path(root, "RNA_Sieve_monaco_pbmcs/standardized/expression.csv"),
    truth = file.path(root, "RNA_Sieve_monaco_pbmcs/standardized/true_proportions.csv"),
    previous_metrics = file.path(
      root,
      "RNA_Sieve_monaco_pbmcs/seven_tools_compare/CIBERSORT_native_ref_top1000_metrics.txt"
    )
  )
)

# LM22 broad-category aggregation. These are direct sums of its 22 fractions.
lm22_map <- list(
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

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

evaluate <- function(pred, truth, method_name) {
  pred <- pred[rownames(truth), , drop = FALSE]
  rows <- lapply(intersect(colnames(truth), names(lm22_map)), function(cell) {
    states <- intersect(lm22_map[[cell]], colnames(pred))
    p <- rowSums(pred[, states, drop = FALSE])
    y <- truth[, cell]
    ok <- is.finite(p) & is.finite(y)
    fit <- if (sum(ok) >= 3 && sd(p[ok]) > 0 && sd(y[ok]) > 0) {
      coef(lm(p[ok] ~ y[ok]))
    } else {
      c(`(Intercept)` = NA_real_, y = NA_real_)
    }
    data.frame(
      Method = method_name,
      CellType = cell,
      PredictedStates = paste(states, collapse = "+"),
      Pearson = safe_cor(p, y, "pearson"),
      Spearman = safe_cor(p, y, "spearman"),
      RMSE = sqrt(mean((p[ok] - y[ok])^2)),
      MAE = mean(abs(p[ok] - y[ok])),
      Slope = unname(fit[2]),
      Intercept = unname(fit[1]),
      N_samples = sum(ok)
    )
  })
  do.call(rbind, rows)
}

summarize_metrics <- function(x, dataset, method) {
  data.frame(
    Dataset = dataset,
    Method = method,
    N_targets = nrow(x),
    ValidPearsonTargets = sum(is.finite(x$Pearson)),
    MeanPearson = mean(x$Pearson, na.rm = TRUE),
    MedianPearson = median(x$Pearson, na.rm = TRUE),
    MeanSpearman = mean(x$Spearman, na.rm = TRUE),
    MeanRMSE = mean(x$RMSE, na.rm = TRUE),
    MeanMAE = mean(x$MAE, na.rm = TRUE)
  )
}

all_metrics <- list()
all_summaries <- list()
for (dataset in names(datasets)) {
  cfg <- datasets[[dataset]]
  bulk <- as.matrix(read.csv(cfg$expression, row.names = 1, check.names = FALSE))
  truth <- as.matrix(read.csv(cfg$truth, row.names = 1, check.names = FALSE))
  storage.mode(bulk) <- "numeric"
  storage.mode(truth) <- "numeric"
  common_samples <- intersect(colnames(bulk), rownames(truth))
  bulk <- bulk[, common_samples, drop = FALSE]
  truth <- truth[common_samples, , drop = FALSE]

  data(LM22, package = "CIBERSORT")
  common_genes <- intersect(rownames(LM22), rownames(bulk))
  message(dataset, ": ", length(common_genes), "/547 LM22 genes; ",
          length(common_samples), " aligned samples")

  # Use the package's default arguments exactly:
  # maxSize=500, perm=0, QN=TRUE.
  fit <- CIBERSORT::cibersort(
    sig_matrix = LM22[common_genes, , drop = FALSE],
    mixture_file = bulk[common_genes, , drop = FALSE]
  )
  fit <- as.matrix(fit)
  fraction_cols <- intersect(colnames(LM22), colnames(fit))
  pred <- fit[, fraction_cols, drop = FALSE]
  pred <- pred[common_samples, , drop = FALSE]

  metrics_lm22 <- evaluate(pred, truth, "CIBERSORT_default_LM22")
  metrics_lm22$Dataset <- dataset
  previous <- read.delim(cfg$previous_metrics, check.names = FALSE)
  previous$Method <- "CIBERSORT_custom_reference_top1000"
  previous$Dataset <- dataset
  previous$N_samples <- previous$N
  previous$N <- NULL
  previous <- previous[, colnames(metrics_lm22), drop = FALSE]

  combined <- rbind(metrics_lm22, previous)
  all_metrics[[dataset]] <- combined
  all_summaries[[paste0(dataset, "_LM22")]] <-
    summarize_metrics(metrics_lm22, dataset, "CIBERSORT_default_LM22")
  all_summaries[[paste0(dataset, "_custom")]] <-
    summarize_metrics(previous, dataset, "CIBERSORT_custom_reference_top1000")

  write.table(pred, file.path(out_dir, paste0(dataset, "_CIBERSORT_default_LM22_fraction.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(metrics_lm22, file.path(out_dir, paste0(dataset, "_CIBERSORT_default_LM22_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(combined, file.path(out_dir, paste0(dataset, "_CIBERSORT_LM22_vs_custom_metrics.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

metrics_all <- do.call(rbind, all_metrics)
summary_all <- do.call(rbind, all_summaries)
rownames(summary_all) <- NULL

wide <- reshape(
  metrics_all[, c("Dataset", "CellType", "Method", "Pearson")],
  idvar = c("Dataset", "CellType"), timevar = "Method", direction = "wide"
)
colnames(wide) <- sub("^Pearson\\.", "", colnames(wide))

write.table(metrics_all, file.path(out_dir, "Newman_Monaco_CIBERSORT_LM22_vs_custom_all_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_all, file.path(out_dir, "Newman_Monaco_CIBERSORT_LM22_vs_custom_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wide, file.path(out_dir, "Newman_Monaco_CIBERSORT_LM22_vs_custom_pearson_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nCIBERSORT default LM22 versus custom-reference CIBERSORT:\n")
print(summary_all, row.names = FALSE)
cat("\nPer-cell-type Pearson:\n")
print(wide, row.names = FALSE)
cat("\nResults written to:", out_dir, "\n")
