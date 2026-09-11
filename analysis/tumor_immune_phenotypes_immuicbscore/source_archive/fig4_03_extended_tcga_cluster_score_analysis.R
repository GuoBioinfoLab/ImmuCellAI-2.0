suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_out <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
out_dir <- file.path(base_out, "extended_TCGA_cluster_score_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(
  fig4_dir,
  "cluster4_class_Unified_ICB_score_association_full",
  "cluster4_class_samples_with_unified_ICB_score.txt"
)
or_file <- file.path(
  fig4_dir,
  "cluster4_class_Unified_ICB_score_association_full",
  "cluster4_class_logistic_OR_per_1SD_score.txt"
)

stopifnot(file.exists(score_file))

cluster_cols <- c(
  "C1" = "#E15759FF",
  "C2" = "#4E79A7FF",
  "C3" = "#F28E2BFF",
  "C4" = "#76B7B2FF"
)

theme_fig4 <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "#F7E7EA", color = "black", linewidth = 0.45),
      strip.text = element_text(color = "black", face = "bold"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", color = "black", hjust = 0),
      plot.subtitle = element_text(color = "grey30", hjust = 0)
    )
}

save_pdf_png <- function(plot, file_stub, width, height) {
  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), plot, width = width, height = height, units = "in")
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), plot, width = width, height = height, units = "in", dpi = 320)
}

auc_manual <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2) return(NA_real_)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  # Positive means values in x are larger than values in y.
  cmp <- outer(x, y, "-")
  (sum(cmp > 0) - sum(cmp < 0)) / (length(x) * length(y))
}

message("Reading TCGA cluster-score table...")
df <- fread(score_file, data.table = FALSE)
df <- df[is.finite(df$Unified_ICB_score) & !is.na(df$cluster), ]
df$cluster <- factor(df$cluster, levels = paste0("C", 1:4))
df$score_z <- as.numeric(scale(df$Unified_ICB_score))

decile_breaks <- quantile(df$Unified_ICB_score, probs = seq(0, 1, 0.1), na.rm = TRUE)
decile_breaks[1] <- decile_breaks[1] - 1e-12
decile_breaks[length(decile_breaks)] <- decile_breaks[length(decile_breaks)] + 1e-12
df$score_decile <- cut(
  df$Unified_ICB_score,
  breaks = decile_breaks,
  labels = paste0("D", 1:10),
  include.lowest = TRUE
)

decile_comp <- as.data.table(df)[, .N, by = .(score_decile, cluster)]
decile_comp[, Total := sum(N), by = score_decile]
decile_comp[, Proportion := N / Total]
decile_comp$score_decile <- factor(decile_comp$score_decile, levels = paste0("D", 1:10))
fwrite(decile_comp, file.path(out_dir, "score_decile_cluster_composition.txt"), sep = "\t")

p_decile <- ggplot(decile_comp, aes(x = score_decile, y = Proportion, fill = cluster)) +
  geom_col(width = 0.82, color = "black", linewidth = 0.18) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "TCGA cluster composition across Unified ICB Score deciles",
    subtitle = "Higher deciles show C4 expansion and C1 contraction",
    x = "Unified ICB Score decile",
    y = "Cluster proportion",
    fill = "Cluster"
  ) +
  theme_fig4(base_size = 10)
save_pdf_png(p_decile, "Figure5_score_decile_cluster_composition", 7.2, 4.2)

tab <- table(df$score_decile, df$cluster)
chisq_res <- suppressWarnings(chisq.test(tab))
residual_df <- as.data.frame(as.table(chisq_res$stdres))
colnames(residual_df) <- c("score_decile", "cluster", "standardized_residual")
expected <- as.data.frame(as.table(chisq_res$expected))
colnames(expected) <- c("score_decile", "cluster", "expected")
observed <- as.data.frame(as.table(tab))
colnames(observed) <- c("score_decile", "cluster", "observed")
residual_df <- merge(residual_df, expected, by = c("score_decile", "cluster"))
residual_df <- merge(residual_df, observed, by = c("score_decile", "cluster"))
residual_df$log2_observed_expected <- log2((residual_df$observed + 0.5) / (residual_df$expected + 0.5))
residual_df$score_decile <- factor(residual_df$score_decile, levels = paste0("D", 1:10))
residual_df$cluster <- factor(residual_df$cluster, levels = paste0("C", 1:4))
fwrite(residual_df, file.path(out_dir, "score_decile_cluster_enrichment_residuals.txt"), sep = "\t")

chisq_summary <- data.frame(
  test = "Cluster composition differs across Unified ICB Score deciles",
  statistic = unname(chisq_res$statistic),
  df = unname(chisq_res$parameter),
  p_value = chisq_res$p.value
)
fwrite(chisq_summary, file.path(out_dir, "score_decile_cluster_composition_chisq.txt"), sep = "\t")

p_heat <- ggplot(residual_df, aes(x = score_decile, y = cluster, fill = standardized_residual)) +
  geom_tile(color = "white", linewidth = 0.55) +
  scale_fill_gradient2(
    low = "#3B74A8",
    mid = "white",
    high = "#B23A48",
    midpoint = 0,
    name = "Std. residual"
  ) +
  labs(
    title = "Cluster enrichment residuals across score deciles",
    subtitle = "Red: observed more than expected; blue: observed less than expected",
    x = "Unified ICB Score decile",
    y = "TCGA immune cluster"
  ) +
  theme_fig4(base_size = 10) +
  theme(panel.grid = element_blank())
save_pdf_png(p_heat, "Figure6_score_decile_cluster_enrichment_residual_heatmap", 6.6, 3.4)

pairwise_rows <- list()
clusters <- levels(df$cluster)
idx <- 1L
for (a in clusters) {
  for (b in clusters) {
    if (which(clusters == a) >= which(clusters == b)) next
    xa <- df$Unified_ICB_score[df$cluster == a]
    xb <- df$Unified_ICB_score[df$cluster == b]
    wt <- suppressWarnings(wilcox.test(xa, xb))
    pairwise_rows[[idx]] <- data.frame(
      comparison = paste0(a, " vs ", b),
      cluster_a = a,
      cluster_b = b,
      mean_a = mean(xa, na.rm = TRUE),
      mean_b = mean(xb, na.rm = TRUE),
      mean_diff_a_minus_b = mean(xa, na.rm = TRUE) - mean(xb, na.rm = TRUE),
      cliffs_delta_a_minus_b = cliffs_delta(xa, xb),
      wilcox_p = wt$p.value,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
pairwise_df <- rbindlist(pairwise_rows)
pairwise_df$BH_p <- p.adjust(pairwise_df$wilcox_p, method = "BH")
fwrite(pairwise_df, file.path(out_dir, "pairwise_cluster_score_effect_size.txt"), sep = "\t")

effect_c4 <- pairwise_df[pairwise_df$cluster_a == "C4" | pairwise_df$cluster_b == "C4", ]
effect_c4$other_cluster <- ifelse(effect_c4$cluster_a == "C4", effect_c4$cluster_b, effect_c4$cluster_a)
effect_c4$C4_minus_other_mean <- ifelse(
  effect_c4$cluster_a == "C4",
  effect_c4$mean_diff_a_minus_b,
  -effect_c4$mean_diff_a_minus_b
)
effect_c4$C4_minus_other_delta <- ifelse(
  effect_c4$cluster_a == "C4",
  effect_c4$cliffs_delta_a_minus_b,
  -effect_c4$cliffs_delta_a_minus_b
)
effect_c4$other_cluster <- factor(effect_c4$other_cluster, levels = c("C1", "C2", "C3"))

p_effect <- ggplot(effect_c4, aes(x = other_cluster, y = C4_minus_other_delta, fill = other_cluster)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_col(width = 0.62, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = cluster_cols[c("C1", "C2", "C3")]) +
  labs(
    title = "Effect size for higher C4 Unified ICB Score",
    subtitle = "Cliff's delta > 0 means C4 has higher score than the compared cluster",
    x = "Compared cluster",
    y = "Cliff's delta: C4 minus other"
  ) +
  theme_fig4(base_size = 10) +
  theme(legend.position = "none")
save_pdf_png(p_effect, "Figure7_C4_vs_other_clusters_score_effect_size", 4.6, 3.7)

if (file.exists(or_file)) {
  or_df <- fread(or_file, data.table = FALSE)
  colnames(or_df) <- sub("^target_cluster$", "cluster", colnames(or_df))
  or_df$cluster <- factor(or_df$cluster, levels = paste0("C", 1:4))
  fwrite(or_df, file.path(out_dir, "logistic_OR_per_1SD_score_for_cluster_membership.txt"), sep = "\t")

  p_or <- ggplot(or_df, aes(x = OR_per_1SD_score, y = cluster, color = cluster)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.18, linewidth = 0.7) +
    geom_point(size = 2.6) +
    scale_color_manual(values = cluster_cols) +
    scale_x_log10() +
    labs(
      title = "Cluster membership odds per 1-SD Unified ICB Score increase",
      subtitle = "OR > 1 indicates enrichment at higher score",
      x = "Odds ratio per 1-SD score increase (log scale)",
      y = "TCGA immune cluster"
    ) +
    theme_fig4(base_size = 10) +
    theme(legend.position = "none")
  save_pdf_png(p_or, "Figure8_cluster_membership_OR_per_1SD_score", 5.6, 3.7)
}

p_ecdf <- ggplot(df, aes(x = Unified_ICB_score, color = cluster)) +
  stat_ecdf(linewidth = 0.85) +
  scale_color_manual(values = cluster_cols) +
  labs(
    title = "Cumulative distribution of Unified ICB Score by cluster",
    subtitle = "Right-shifted curves indicate higher score distribution",
    x = "Unified ICB Score",
    y = "Empirical cumulative probability",
    color = "Cluster"
  ) +
  theme_fig4(base_size = 10)
save_pdf_png(p_ecdf, "Figure9_cluster_score_ECDF", 5.5, 4.0)

cluster_auc_rows <- list()
for (cl in levels(df$cluster)) {
  cluster_auc_rows[[cl]] <- data.frame(
    cluster = cl,
    AUC_high_score_predicts_cluster = auc_manual(as.integer(df$cluster == cl), df$Unified_ICB_score),
    N_cluster = sum(df$cluster == cl),
    N_total = nrow(df),
    stringsAsFactors = FALSE
  )
}
cluster_auc <- rbindlist(cluster_auc_rows)
cluster_auc$Direction <- ifelse(
  cluster_auc$AUC_high_score_predicts_cluster >= 0.5,
  "higher score enriched",
  "lower score enriched"
)
cluster_auc$OrientedAUC <- pmax(cluster_auc$AUC_high_score_predicts_cluster, 1 - cluster_auc$AUC_high_score_predicts_cluster)
fwrite(cluster_auc, file.path(out_dir, "score_predict_cluster_membership_auc.txt"), sep = "\t")

run_info <- data.frame(
  item = c("input_score_file", "n_samples", "chisq_p_decile_cluster", "analysis_date"),
  value = c(score_file, nrow(df), signif(chisq_res$p.value, 4), as.character(Sys.Date()))
)
fwrite(run_info, file.path(out_dir, "extended_TCGA_cluster_score_analysis_run_info.txt"), sep = "\t")

pdf(file.path(out_dir, "Figure_extended_TCGA_cluster_score_analysis_multipage.pdf"), width = 7.2, height = 4.4)
print(p_decile)
print(p_heat)
print(p_effect)
if (exists("p_or")) print(p_or)
print(p_ecdf)
dev.off()

message("Extended TCGA cluster-score analyses written to: ", out_dir)
