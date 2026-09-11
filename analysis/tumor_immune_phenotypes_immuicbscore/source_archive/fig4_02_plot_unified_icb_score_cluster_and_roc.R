suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
out_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_score_dir <- file.path(fig4_dir, "cluster4_class_Unified_ICB_score_association_full")
model_dir <- file.path(
  fig4_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "more_unified_ICB_score_models"
)

score_sample_file <- file.path(cluster_score_dir, "cluster4_class_samples_with_unified_ICB_score.txt")
quartile_file <- file.path(cluster_score_dir, "cluster4_class_composition_by_score_quartile.txt")
delta_file <- file.path(cluster_score_dir, "cluster4_class_top_vs_bottom_score_enrichment.txt")
auc_file <- file.path(model_dir, "BEST_by_LOSO_auc_by_study.txt")
prediction_file <- file.path(model_dir, "LOSO_predictions.txt")

stopifnot(file.exists(score_sample_file), file.exists(quartile_file), file.exists(delta_file))
stopifnot(file.exists(auc_file), file.exists(prediction_file))

cluster_cols <- c(
  "C1" = "#E15759FF",
  "C2" = "#4E79A7FF",
  "C3" = "#F28E2BFF",
  "C4" = "#76B7B2FF"
)

theme_fig4 <- function(base_size = 11) {
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

roc_coordinates <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2) return(NULL)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  score <- score[ord]
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  data.frame(
    threshold = c(Inf, score, -Inf),
    FPR = c(0, cumsum(y == 0) / n_neg, 1),
    TPR = c(0, cumsum(y == 1) / n_pos, 1),
    stringsAsFactors = FALSE
  )
}

message("Reading TCGA cluster-score association data...")
score_df <- fread(score_sample_file, data.table = FALSE)
quartile_df <- fread(quartile_file, data.table = FALSE)
delta_df <- fread(delta_file, data.table = FALSE)

score_df$cluster <- factor(score_df$cluster, levels = paste0("C", 1:4))
quartile_df$cluster <- factor(quartile_df$cluster, levels = paste0("C", 1:4))
quartile_df$score_quartile <- factor(
  quartile_df$score_quartile,
  levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
)
delta_df$target_cluster <- factor(delta_df$target_cluster, levels = paste0("C", 1:4))

cluster_summary <- as.data.table(score_df)[, .(
  N = .N,
  MeanScore = mean(Unified_ICB_score, na.rm = TRUE),
  MedianScore = median(Unified_ICB_score, na.rm = TRUE)
), by = cluster][order(cluster)]
fwrite(cluster_summary, file.path(out_dir, "TCGA_cluster_Unified_ICB_score_summary.txt"), sep = "\t")

p_score <- ggplot(score_df, aes(x = cluster, y = Unified_ICB_score, fill = cluster)) +
  geom_violin(width = 0.86, color = "black", linewidth = 0.35, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.size = 0.35, color = "black", linewidth = 0.35, alpha = 0.88) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2, fill = "white", color = "black") +
  scale_fill_manual(values = cluster_cols) +
  labs(
    title = "Unified ICB Score across TCGA immune clusters",
    subtitle = "Higher values indicate a higher predicted probability of ICB response",
    x = "TCGA immune cluster",
    y = "Unified ICB Score"
  ) +
  theme_fig4() +
  theme(legend.position = "none")
save_pdf_png(p_score, "Figure1_TCGA_cluster_Unified_ICB_score_violin", 4.8, 4.0)

target_quartile <- quartile_df[quartile_df$cluster %in% c("C1", "C3", "C4"), ]
p_quartile <- ggplot(target_quartile, aes(x = score_quartile, y = Proportion, group = cluster, color = cluster)) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.4) +
  scale_color_manual(values = cluster_cols[c("C1", "C3", "C4")]) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Cluster composition across score quartiles",
    subtitle = "C4 enrichment and C1 depletion from low to high Unified ICB Score",
    x = "Unified ICB Score quartile",
    y = "Cluster proportion",
    color = "Cluster"
  ) +
  theme_fig4() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_pdf_png(p_quartile, "Figure2_C1_C3_C4_proportion_across_score_quartiles", 6.8, 4.1)

p_delta <- ggplot(delta_df, aes(x = target_cluster, y = Delta_Q4_minus_Q1, fill = target_cluster)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_col(width = 0.65, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "High-score enrichment of TCGA immune clusters",
    subtitle = "Positive values indicate enrichment in Q4 highest score versus Q1 lowest score",
    x = "TCGA immune cluster",
    y = "Q4 - Q1 cluster proportion"
  ) +
  theme_fig4() +
  theme(legend.position = "none")
save_pdf_png(p_delta, "Figure3_Q4_vs_Q1_cluster_enrichment", 4.8, 4.0)

message("Reading external ICB model predictions...")
auc_df <- fread(auc_file, data.table = FALSE)
pred_df <- fread(prediction_file, data.table = FALSE)

selected_model <- "cell_clr53__ranger_balanced"
high_auc <- auc_df[
  auc_df$model == selected_model & is.finite(auc_df$AUC) & auc_df$AUC >= 0.75,
]
high_auc <- high_auc[order(-high_auc$AUC), ]
fwrite(high_auc, file.path(out_dir, "selected_high_prediction_ICB_cohorts_AUC.txt"), sep = "\t")

pred_sel <- pred_df[
  pred_df$model == selected_model & pred_df$SRA_study %in% high_auc$SRA_study,
]

roc_list <- list()
for (study in unique(pred_sel$SRA_study)) {
  dd <- pred_sel[pred_sel$SRA_study == study, ]
  rr <- roc_coordinates(dd$y, dd$prediction)
  if (!is.null(rr)) {
    rr$SRA_study <- study
    rr$AUC <- auc_manual(dd$y, dd$prediction)
    rr$N <- nrow(dd)
    rr$N_R <- sum(dd$y == 1)
    rr$N_NR <- sum(dd$y == 0)
    roc_list[[study]] <- rr
  }
}
roc_df <- rbindlist(roc_list, fill = TRUE)
roc_df$label <- sprintf("%s\nAUC = %.3f", roc_df$SRA_study, roc_df$AUC)
label_order <- unique(roc_df[order(-roc_df$AUC), c("label", "AUC")])$label
roc_df$label <- factor(roc_df$label, levels = label_order)
fwrite(roc_df, file.path(out_dir, "ROC_coordinates_selected_high_prediction_ICB_cohorts.txt"), sep = "\t")

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_path(color = "#B23A48", linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~ label, ncol = 3) +
  labs(
    title = "External ICB cohorts with strong prediction",
    subtitle = "Selected cohorts: LOSO AUC >= 0.75 for cell_clr53__ranger_balanced",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(base_size = 9) +
  theme(
    strip.text = element_text(color = "black", face = "bold", size = 8.5),
    axis.text = element_text(color = "black", size = 8)
  )
save_pdf_png(p_roc, "Figure4_selected_external_ICB_cohort_ROC", 8.6, 5.6)

pdf(file.path(out_dir, "Figure_Unified_ICB_score_cluster_ICB_integration_multipage.pdf"), width = 8, height = 5.4)
print(p_score)
print(p_quartile)
print(p_delta)
print(p_roc)
dev.off()

run_info <- data.frame(
  item = c(
    "score_model",
    "score_features",
    "tcga_cluster_file_source",
    "external_icb_auc_threshold",
    "n_high_auc_cohorts"
  ),
  value = c(
    selected_model,
    "CLR-transformed 53 ImmuCellAI2 cell fractions",
    score_sample_file,
    "LOSO AUC >= 0.75",
    nrow(high_auc)
  )
)
fwrite(run_info, file.path(out_dir, "Fig4_Unified_ICB_score_integration_run_info.txt"), sep = "\t")

message("Fig4 Unified ICB Score integration figures written to: ", out_dir)
