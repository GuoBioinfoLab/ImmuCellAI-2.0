# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
out_dir <- file.path(work_dir, "Fig2", "mean_pearson_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- data.frame(
  Dataset = c("GSE146771", "GSE164522", "GSE176078", "GSE107011"),
  File = c(
    fig3_tumor_metrics("GSE146771"),
    fig3_tumor_metrics("GSE164522"),
    fig3_tumor_metrics("GSE176078"),
    file.path(work_dir, "GSE107011_adjusted_7tools_flatvb/GSE107011_adjusted_7tools_per_celltype_metrics.txt")
  ),
  stringsAsFactors = FALSE
)

method_labels <- c(
  ImmuCellAI2_tcell_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_flat_VB_UNKNOWN = "ImmuCellAI 2.0",
  ImmuCellAI2_tcell_VB_UNKNOWN = "ImmuCellAI 2.0",
  ImmuCellAI2_best_ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_best_ImmuCellAI2_flat_VB_UNKNOWN = "ImmuCellAI 2.0",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  DWLS_weighted_lm = "DWLS",
  MuSiC_basic = "MuSiC",
  CIBERSORT_default_LM22 = "CIBERSORT",
  CITMIC_native = "CITMIC",
  CIBERSORT_nuSVR_ref_top1000 = "CIBERSORT",
  ImmuCellAI_native = "ImmuCellAI"
)

method_order <- c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC")
dataset_order <- c("GSE146771", "GSE164522", "GSE176078", "GSE107011")

read_one_metrics <- function(dataset, file) {
  if (!file.exists(file)) stop("Input file not found: ", file)
  x <- read.delim(file, check.names = FALSE)
  required <- c("Method", "CellType", "Pearson", "Spearman", "RMSE", "MAE")
  miss <- setdiff(required, colnames(x))
  if (length(miss) > 0) stop("Missing columns in ", file, ": ", paste(miss, collapse = ", "))
  x$Dataset <- dataset
  x$Method <- as.character(x$Method)
  x$CellType <- as.character(x$CellType)
  x$Pearson <- suppressWarnings(as.numeric(x$Pearson))
  x$Spearman <- suppressWarnings(as.numeric(x$Spearman))
  x$RMSE <- suppressWarnings(as.numeric(x$RMSE))
  x$MAE <- suppressWarnings(as.numeric(x$MAE))
  x$Slope <- if ("Slope" %in% colnames(x)) suppressWarnings(as.numeric(x$Slope)) else NA_real_
  x$Intercept <- if ("Intercept" %in% colnames(x)) suppressWarnings(as.numeric(x$Intercept)) else NA_real_
  x$MethodShort <- ifelse(x$Method %in% names(method_labels), unname(method_labels[x$Method]), x$Method)
  keep_cols <- c("Dataset", "Method", "MethodShort", "CellType", "Pearson", "Spearman", "RMSE", "MAE", "Slope", "Intercept")
  x[, keep_cols, drop = FALSE]
}

all_metrics <- do.call(
  rbind,
  Map(read_one_metrics, input_files$Dataset, input_files$File)
)

# Match the current pie-plot filtering rules.
drop_rules <- data.frame(
  Dataset = c(rep("GSE176078", 3), "GSE164522"),
  CellType = c("NKT", "Th1", "Th17", "NKT"),
  stringsAsFactors = FALSE
)
drop_key <- paste(drop_rules$Dataset, drop_rules$CellType)
all_metrics <- all_metrics[!paste(all_metrics$Dataset, all_metrics$CellType) %in% drop_key, , drop = FALSE]
all_metrics <- all_metrics[all_metrics$MethodShort %in% method_order, , drop = FALSE]

summarize_one <- function(x) {
  data.frame(
    N_targets = length(unique(x$CellType)),
    MeanPearson = mean(x$Pearson, na.rm = TRUE),
    MedianPearson = median(x$Pearson, na.rm = TRUE),
    MeanSpearman = mean(x$Spearman, na.rm = TRUE),
    MedianSpearman = median(x$Spearman, na.rm = TRUE),
    MeanRMSE = mean(x$RMSE, na.rm = TRUE),
    MedianRMSE = median(x$RMSE, na.rm = TRUE),
    MeanMAE = mean(x$MAE, na.rm = TRUE),
    MedianMAE = median(x$MAE, na.rm = TRUE),
    MeanSlope = mean(x$Slope, na.rm = TRUE),
    MedianSlope = median(x$Slope, na.rm = TRUE),
    MeanSlopeAbsError = abs(mean(x$Slope, na.rm = TRUE) - 1),
    MedianSlopeAbsError = abs(median(x$Slope, na.rm = TRUE) - 1),
    MeanIntercept = mean(x$Intercept, na.rm = TRUE),
    MedianIntercept = median(x$Intercept, na.rm = TRUE),
    MeanAbsIntercept = abs(mean(x$Intercept, na.rm = TRUE)),
    MedianAbsIntercept = abs(median(x$Intercept, na.rm = TRUE)),
    ValidPearsonTargets = sum(!is.na(x$Pearson)),
    Coverage = sum(!is.na(x$Pearson)) / length(unique(x$CellType)),
    PenalizedMeanPearson = mean(ifelse(is.na(x$Pearson), -1, x$Pearson)),
    stringsAsFactors = FALSE
  )
}

keys <- unique(all_metrics[, c("Dataset", "MethodShort")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  sub <- all_metrics[all_metrics$Dataset == keys$Dataset[i] & all_metrics$MethodShort == keys$MethodShort[i], , drop = FALSE]
  cbind(Dataset = keys$Dataset[i], Method = keys$MethodShort[i], summarize_one(sub))
}))
summary_table$Dataset <- factor(summary_table$Dataset, levels = dataset_order)
summary_table$Method <- factor(summary_table$Method, levels = method_order)
summary_table <- summary_table[order(summary_table$Dataset, summary_table$Method), ]

write.table(all_metrics, file.path(out_dir, "all_metrics_current_plot_scope.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_table, file.path(out_dir, "overall_summary_current_plot_scope.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package RColorBrewer is required.")
library(ggplot2)

method_colors <- c(
  "ImmuCellAI 2.0" = "#F66463",
  BayesPrism = "#379DA5",
  DWLS = "#6BA5C3",
  MuSiC = "#FAC74C",
  CITMIC = "#95A8AC",
  CIBERSORT = "#AE997E",
  ImmuCellAI = "#8E6BBE"
)

theme_blue <- theme(
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "gray80", linewidth = 0.6, linetype = "dashed"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.8),
  strip.text = element_text(size = 13, colour = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 11, colour = "black"),
  axis.text.y = element_text(size = 11, colour = "black"),
  axis.title = element_text(size = 12, face = "bold", colour = "darkred"),
  legend.title = element_text(size = 11, colour = "black"),
  legend.text = element_text(size = 11, colour = "black"),
  legend.key = element_blank()
)

to_long <- function(df, metrics) {
  out <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(
      Dataset = df$Dataset,
      Method = df$Method,
      Metric = m,
      Value = suppressWarnings(as.numeric(df[[m]])),
      stringsAsFactors = FALSE
    )
  }))
  out$Dataset <- factor(out$Dataset, levels = dataset_order)
  out$Method <- factor(out$Method, levels = method_order)
  out
}

higher_metrics <- c("MeanPearson")
lower_metrics <- c("MeanRMSE", "MeanMAE")

higher_long <- to_long(summary_table, higher_metrics)
lower_long <- to_long(summary_table, lower_metrics)

p_higher <- ggplot(higher_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Higher is better") +
  theme_blue +
  theme(legend.position = "bottom")

p_lower <- ggplot(lower_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Lower is better") +
  theme_blue +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "overall_mean_pearson_bar.pdf"), p_higher, width = 18, height = 4.8, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_mean_pearson_bar.png"), p_higher, width = 18, height = 4.8, dpi = 220, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_lower_better_error_bar.pdf"), p_lower, width = 18, height = 6, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_lower_better_error_bar.png"), p_lower, width = 18, height = 6, dpi = 220, limitsize = FALSE)

# A compact single-metric view for the main text.
p_penalized <- ggplot(summary_table, aes(x = Method, y = PenalizedMeanPearson, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(. ~ Dataset) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Penalized mean Pearson") +
  theme_blue +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "overall_penalized_mean_pearson_bar.pdf"), p_penalized, width = 18, height = 4.8, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_penalized_mean_pearson_bar.png"), p_penalized, width = 18, height = 4.8, dpi = 220, limitsize = FALSE)

slope_long <- to_long(summary_table, c("MeanSlope", "MedianSlope"))
p_slope <- ggplot(slope_long, aes(x = Method, y = Value, fill = Method)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Slope; 1 is ideal") +
  theme_blue +
  theme(legend.position = "bottom")

slope_error_long <- to_long(summary_table, c("MeanSlopeAbsError", "MedianSlopeAbsError"))
p_slope_error <- ggplot(slope_error_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Absolute distance from slope = 1; lower is better") +
  theme_blue +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "overall_slope_bar.pdf"), p_slope, width = 18, height = 6, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_slope_bar.png"), p_slope, width = 18, height = 6, dpi = 220, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_slope_abs_error_bar.pdf"), p_slope_error, width = 18, height = 6, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_slope_abs_error_bar.png"), p_slope_error, width = 18, height = 6, dpi = 220, limitsize = FALSE)

intercept_long <- to_long(summary_table, c("MeanIntercept", "MedianIntercept"))
p_intercept <- ggplot(intercept_long, aes(x = Method, y = Value, fill = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Intercept; 0 is ideal") +
  theme_blue +
  theme(legend.position = "bottom")

intercept_error_long <- to_long(summary_table, c("MeanAbsIntercept", "MedianAbsIntercept"))
p_intercept_error <- ggplot(intercept_error_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  facet_grid(Metric ~ Dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = "Absolute intercept; lower is better") +
  theme_blue +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "overall_intercept_bar.pdf"), p_intercept, width = 18, height = 6, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_intercept_bar.png"), p_intercept, width = 18, height = 6, dpi = 220, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_abs_intercept_bar.pdf"), p_intercept_error, width = 18, height = 6, limitsize = FALSE)
ggsave(file.path(out_dir, "overall_abs_intercept_bar.png"), p_intercept_error, width = 18, height = 6, dpi = 220, limitsize = FALSE)

print(summary_table)
message("Saved output directory: ", out_dir)
