options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
score_cluster_dir <- file.path(fig4_dir, "Unified_ICB_score_TCGA_cluster_association")
out_dir <- file.path(fig4_dir, "Unified_ICB_score_C3_C4_axis_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcga_fraction_file <- file.path(
  fig4_dir,
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000",
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
tcga_cluster_file <- file.path(cluster_dir, "cluster4_class.txt")
tcga_score_file <- file.path(score_cluster_dir, "TCGA_samples_unified_ICB_score_by_cluster.txt")
icb_axis_file <- file.path(base_dir, "C3_C4_axis_ICB_response_scores", "external_ICB_C3_C4_axis_scores.txt")
icb_score_file <- file.path(base_dir, "more_unified_ICB_score_models", "LOSO_predictions.txt")
icb_meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(data.frame(correlation = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(correlation = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
}

read_fraction <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  rownames(x) <- as.character(x[[1]])
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m
}

clr_transform <- function(x, eps = 1e-5) {
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx), "-")
}

message("Reading TCGA score and rebuilding C3/C4 axis from cluster4_class.txt...")
tcga_fraction <- read_fraction(tcga_fraction_file)
tcga_cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
tcga_cluster <- tcga_cluster_raw[, c("sample", "cluster")]
tcga_cluster$sample <- as.character(tcga_cluster$sample)
tcga_cluster$cluster <- paste0("C", as.integer(tcga_cluster$cluster))
tcga_cluster <- tcga_cluster[tcga_cluster$sample %in% rownames(tcga_fraction), , drop = FALSE]
tcga_fraction <- tcga_fraction[tcga_cluster$sample, , drop = FALSE]

tcga_score <- fread(tcga_score_file, data.table = FALSE, check.names = FALSE)
tcga_score <- tcga_score[tcga_score$sample %in% rownames(tcga_fraction), , drop = FALSE]
tcga_fraction <- tcga_fraction[tcga_score$sample, , drop = FALSE]

tcga_clr <- clr_transform(tcga_fraction)
centroids_raw <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- intersect(tcga_cluster$sample[tcga_cluster$cluster == cl], rownames(tcga_fraction))
  colMeans(tcga_fraction[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_raw) <- paste0("C", 1:4)

centroids_clr <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- intersect(tcga_cluster$sample[tcga_cluster$cluster == cl], rownames(tcga_clr))
  colMeans(tcga_clr[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_clr) <- paste0("C", 1:4)

axis_tcga <- data.frame(sample = rownames(tcga_fraction), stringsAsFactors = FALSE)
for (cl in paste0("C", 1:4)) {
  axis_tcga[[paste0(cl, "_spearman_raw")]] <- apply(tcga_fraction, 1, function(v) {
    suppressWarnings(cor(v, centroids_raw[cl, ], method = "spearman"))
  })
  axis_tcga[[paste0(cl, "_pearson_clr")]] <- apply(tcga_clr, 1, function(v) {
    suppressWarnings(cor(v, centroids_clr[cl, ], method = "pearson"))
  })
}
axis_tcga$C3_low_clr <- -axis_tcga$C3_pearson_clr
axis_tcga$C4_minus_C3_spearman_raw <- axis_tcga$C4_spearman_raw - axis_tcga$C3_spearman_raw
axis_tcga$C4_minus_C3_pearson_clr <- axis_tcga$C4_pearson_clr - axis_tcga$C3_pearson_clr
axis_tcga$C4_minus_C2_spearman_raw <- axis_tcga$C4_spearman_raw - axis_tcga$C2_spearman_raw

tcga_df <- merge(tcga_score, axis_tcga, by = "sample", all.x = TRUE)
tcga_df$Dataset <- "TCGA"

message("Reading ICB LOSO score and C3/C4 axis scores...")
icb_axis <- fread(icb_axis_file, data.table = FALSE, check.names = FALSE)
icb_pred <- fread(icb_score_file, data.table = FALSE, check.names = FALSE)
icb_pred <- icb_pred[icb_pred$model == "cell_clr53__ranger_balanced", , drop = FALSE]
icb_meta <- fread(icb_meta_file, data.table = FALSE, check.names = FALSE)

response <- rep(NA_integer_, nrow(icb_meta))
response[icb_meta$ResponseGroup == "Responder"] <- 1L
response[icb_meta$ResponseGroup == "Non-responder"] <- 0L
response[is.na(response) & icb_meta$Response == "R"] <- 1L
response[is.na(response) & icb_meta$Response == "NR"] <- 0L
icb_meta$response_binary <- response
icb_meta <- icb_meta[!is.na(icb_meta$response_binary), , drop = FALSE]

icb_df <- merge(
  icb_pred[, c("sample", "SRA_study", "y", "prediction")],
  icb_axis,
  by = "sample",
  all.x = TRUE
)
icb_df <- merge(
  icb_df,
  icb_meta[, c("Run", "Cancer", "Cancer_type", "Drug", "disease", "response_binary")],
  by.x = "sample",
  by.y = "Run",
  all.x = TRUE
)
icb_df$Unified_ICB_score <- icb_df$prediction
icb_df$C3_low_clr <- -icb_df$C3_pearson_clr
icb_df$Dataset <- "ICB"

axis_vars <- c(
  "C1_spearman_raw", "C2_spearman_raw", "C3_spearman_raw", "C4_spearman_raw",
  "C1_pearson_clr", "C2_pearson_clr", "C3_pearson_clr", "C4_pearson_clr",
  "C3_low_clr",
  "C4_minus_C3_spearman_raw", "C4_minus_C3_pearson_clr",
  "C4_minus_C2_spearman_raw"
)
axis_vars <- unique(axis_vars)

make_cor_table <- function(df, score_col, dataset_name) {
  rows <- list()
  for (v in intersect(axis_vars, colnames(df))) {
    for (method in c("spearman", "pearson")) {
      tmp <- safe_cor_test(df[[score_col]], df[[v]], method = method)
      rows[[length(rows) + 1L]] <- data.frame(
        Dataset = dataset_name,
        score = score_col,
        axis_variable = v,
        method = method,
        tmp,
        stringsAsFactors = FALSE
      )
    }
  }
  rbindlist(rows)
}

tcga_cor <- make_cor_table(tcga_df, "Unified_ICB_score", "TCGA")
icb_cor <- make_cor_table(icb_df, "Unified_ICB_score", "ICB")
cor_all <- rbindlist(list(tcga_cor, icb_cor), fill = TRUE)
cor_all$BH_P <- p.adjust(cor_all$p_value, method = "BH")
cor_all <- cor_all[order(cor_all$Dataset, cor_all$method, -abs(cor_all$correlation)), ]

fwrite(tcga_df, file.path(out_dir, "TCGA_unified_score_C3_C4_axis_values.txt"), sep = "\t")
fwrite(icb_df, file.path(out_dir, "ICB_unified_score_C3_C4_axis_values.txt"), sep = "\t")
fwrite(cor_all, file.path(out_dir, "Unified_ICB_score_C3_C4_axis_correlation_table.txt"), sep = "\t")

plot_cor <- cor_all[cor_all$method == "spearman", ]
plot_cor$axis_variable <- factor(plot_cor$axis_variable, levels = rev(unique(plot_cor$axis_variable)))
p1 <- ggplot(plot_cor, aes(x = Dataset, y = axis_variable, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  scale_fill_gradient2(low = "#3341A3", mid = "white", high = "#E15759", midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = NULL, fill = "Spearman r",
       title = "Correlation between Unified ICB score and C3/C4-axis features") +
  theme(panel.grid = element_blank())
ggsave(file.path(out_dir, "Unified_ICB_score_C3_C4_axis_correlation_heatmap.pdf"), p1, width = 6.5, height = 5)
ggsave(file.path(out_dir, "Unified_ICB_score_C3_C4_axis_correlation_heatmap.png"), p1, width = 6.5, height = 5, dpi = 300)

scatter_specs <- data.frame(
  xvar = c("C3_low_clr", "C4_minus_C3_spearman_raw", "C4_pearson_clr", "C3_pearson_clr"),
  label = c("C3-low axis", "C4-C3 Spearman axis", "C4-like CLR", "C3-like CLR"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(scatter_specs))) {
  xvar <- scatter_specs$xvar[i]
  dd <- rbindlist(list(
    data.frame(Dataset = "TCGA", Unified_ICB_score = tcga_df$Unified_ICB_score, axis_value = tcga_df[[xvar]]),
    data.frame(Dataset = "ICB", Unified_ICB_score = icb_df$Unified_ICB_score, axis_value = icb_df[[xvar]])
  ))
  p <- ggplot(dd, aes(x = axis_value, y = Unified_ICB_score)) +
    geom_point(size = 0.55, alpha = 0.45, color = "#379DA5") +
    geom_smooth(method = "lm", se = TRUE, color = "#F66463", linewidth = 0.55) +
    facet_wrap(~ Dataset, scales = "free") +
    theme_bw(base_size = 10) +
    labs(x = scatter_specs$label[i], y = "Unified ICB score",
         title = paste0("Unified ICB score vs ", scatter_specs$label[i])) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, paste0("scatter_Unified_ICB_score_vs_", xvar, ".pdf")), p, width = 7.2, height = 3.6)
  ggsave(file.path(out_dir, paste0("scatter_Unified_ICB_score_vs_", xvar, ".png")), p, width = 7.2, height = 3.6, dpi = 300)
}

top_lines <- capture.output(print(cor_all[cor_all$method == "spearman", ][order(Dataset, -abs(correlation)), ][, c("Dataset", "axis_variable", "correlation", "p_value", "BH_P", "n")], row.names = FALSE))
writeLines(c(
  "Unified ICB score vs C3/C4-axis correlation analysis",
  "TCGA score: ranger ICB-response model projected to TCGA samples.",
  "ICB score: leave-one-study-out RF-CLR prediction score.",
  "",
  top_lines
), file.path(out_dir, "README_Unified_ICB_score_C3_C4_axis_correlation.txt"))

print(cor_all[cor_all$method == "spearman", ][order(Dataset, -abs(correlation)), ])
