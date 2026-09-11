suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
out_dir <- file.path(base_dir, "fivefold_cv_unified_icb_score")

pred_file <- file.path(out_dir, "test_add_high_within_cv_cohorts_expanded_predictions.txt")
auc_file <- file.path(out_dir, "test_add_high_within_cv_cohorts_expanded_per_study_auc.txt")

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
  data.table(
    threshold = c(Inf, score, -Inf),
    FPR = c(0, cumsum(y == 0) / n_neg, 1),
    TPR = c(0, cumsum(y == 1) / n_pos, 1)
  )
}

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

pred <- fread(pred_file, data.table = TRUE)
pred[, ResponseBinary := as.integer(ResponseBinary)]

pooled_auc <- auc_manual(pred$ResponseBinary, pred$Unified_ICB_Score)
pooled_roc <- roc_coordinates(pred$ResponseBinary, pred$Unified_ICB_Score)
pooled_roc[, Curve := sprintf("Pooled 7 cohorts (AUC = %.3f)", pooled_auc)]
fwrite(pooled_roc, file.path(out_dir, "expanded_7cohort_pooled_5fold_ROC_coordinates.txt"), sep = "\t")

label_text <- sprintf(
  "AUC = %.3f\nN = %d\nR = %d, NR = %d",
  pooled_auc,
  nrow(pred),
  sum(pred$ResponseBinary == 1),
  sum(pred$ResponseBinary == 0)
)

p_pooled <- ggplot(pooled_roc, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_path(color = "#B2182B", linewidth = 1.05) +
  annotate("text", x = 0.62, y = 0.18, label = label_text, hjust = 0, size = 4.2, fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Unified ICB Score ROC",
    subtitle = "Expanded 7 selected ICB cohorts, pooled five-fold cross-validation",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(11)

ggsave(file.path(out_dir, "Figure_expanded_7cohort_pooled_5fold_ROC.pdf"), p_pooled, width = 4.8, height = 4.5)
ggsave(file.path(out_dir, "Figure_expanded_7cohort_pooled_5fold_ROC.png"), p_pooled, width = 4.8, height = 4.5, dpi = 500)

per_study_auc <- fread(auc_file, data.table = TRUE)
study_rocs <- rbindlist(lapply(split(pred, pred$SRA_study), function(dd) {
  rr <- roc_coordinates(dd$ResponseBinary, dd$Unified_ICB_Score)
  if (is.null(rr)) return(NULL)
  aa <- auc_manual(dd$ResponseBinary, dd$Unified_ICB_Score)
  rr[, SRA_study := unique(dd$SRA_study)]
  rr[, Label := sprintf("%s\nAUC = %.3f, N = %d", unique(dd$SRA_study), aa, nrow(dd))]
  rr
}), fill = TRUE)

label_levels <- per_study_auc[order(-AUC), SRA_study]
label_map <- unique(study_rocs[, .(SRA_study, Label)])
label_map <- label_map[match(label_levels, SRA_study)]
study_rocs[, Label := factor(Label, levels = label_map$Label)]
fwrite(study_rocs, file.path(out_dir, "expanded_7cohort_per_study_5fold_ROC_coordinates.txt"), sep = "\t")

p_study <- ggplot(study_rocs, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_path(color = "#379DA5", linewidth = 0.75) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~ Label, ncol = 4) +
  labs(
    title = "Per-cohort ROC curves",
    subtitle = "Expanded 7 selected ICB cohorts, pooled five-fold cross-validation predictions",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(9)

ggsave(file.path(out_dir, "Figure_expanded_7cohort_per_study_5fold_ROC.pdf"), p_study, width = 9.8, height = 5.6)
ggsave(file.path(out_dir, "Figure_expanded_7cohort_per_study_5fold_ROC.png"), p_study, width = 9.8, height = 5.6, dpi = 500)

message("Pooled AUC: ", sprintf("%.3f", pooled_auc))
message("ROC figures written to: ", out_dir)
