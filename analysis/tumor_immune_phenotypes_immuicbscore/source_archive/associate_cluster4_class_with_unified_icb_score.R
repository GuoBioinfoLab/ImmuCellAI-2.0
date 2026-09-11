options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
score_dir <- file.path(fig4_dir, "Unified_ICB_score_TCGA_cluster_association")
out_dir <- file.path(fig4_dir, "cluster4_class_Unified_ICB_score_association")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(score_dir, "TCGA_samples_unified_ICB_score_by_cluster.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

auc_manual <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

or_from_glm <- function(df, target_cluster) {
  dat <- data.frame(
    y = as.integer(df$cluster == target_cluster),
    score_z = as.numeric(scale(df$Unified_ICB_score))
  )
  fit <- glm(y ~ score_z, data = dat, family = binomial())
  sm <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  data.frame(
    target_cluster = target_cluster,
    OR_per_1SD_score = exp(coef(fit)["score_z"]),
    Lower95 = exp(ci["score_z", 1]),
    Upper95 = exp(ci["score_z", 2]),
    P = sm["score_z", "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

fisher_top_bottom <- function(df, target_cluster) {
  dd <- df[df$score_quartile %in% c("Q1_lowest", "Q4_highest"), ]
  tab <- table(
    HighScore = dd$score_quartile == "Q4_highest",
    TargetCluster = dd$cluster == target_cluster
  )
  ft <- fisher.test(tab)
  data.frame(
    target_cluster = target_cluster,
    comparison = "Q4_highest_vs_Q1_lowest",
    OR = unname(ft$estimate),
    P = ft$p.value,
    Q1_target_n = sum(dd$score_quartile == "Q1_lowest" & dd$cluster == target_cluster),
    Q1_total_n = sum(dd$score_quartile == "Q1_lowest"),
    Q4_target_n = sum(dd$score_quartile == "Q4_highest" & dd$cluster == target_cluster),
    Q4_total_n = sum(dd$score_quartile == "Q4_highest"),
    stringsAsFactors = FALSE
  )
}

message("Reading Unified ICB score projected to TCGA cluster4_class samples...")
df <- fread(score_file, data.table = FALSE, check.names = FALSE)
df$cluster <- factor(as.character(df$cluster), levels = paste0("C", 1:4))

if (!all(c("sample", "cluster", "Unified_ICB_score") %in% colnames(df))) {
  stop("Input file must contain sample, cluster, and Unified_ICB_score columns.")
}

df <- df[is.finite(df$Unified_ICB_score) & !is.na(df$cluster), , drop = FALSE]

q <- quantile(df$Unified_ICB_score, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
q[1] <- q[1] - 1e-12
q[5] <- q[5] + 1e-12
df$score_quartile <- cut(
  df$Unified_ICB_score,
  breaks = q,
  labels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
  include.lowest = TRUE
)

cluster_summary <- as.data.table(df)[, .(
  N = .N,
  MeanScore = mean(Unified_ICB_score),
  MedianScore = median(Unified_ICB_score),
  SDScore = sd(Unified_ICB_score),
  Q25 = quantile(Unified_ICB_score, 0.25),
  Q75 = quantile(Unified_ICB_score, 0.75)
), by = cluster][order(cluster)]

quartile_counts <- as.data.table(df)[, .N, by = .(score_quartile, cluster)]
quartile_counts[, Total := sum(N), by = score_quartile]
quartile_counts[, Proportion := N / Total]
quartile_counts <- quartile_counts[order(score_quartile, cluster)]

cluster_counts <- table(df$score_quartile, df$cluster)
chisq_cluster_quartile <- suppressWarnings(chisq.test(cluster_counts))
chisq_df <- data.frame(
  test = "Chi-square: cluster composition across score quartiles",
  statistic = unname(chisq_cluster_quartile$statistic),
  df = unname(chisq_cluster_quartile$parameter),
  p_value = chisq_cluster_quartile$p.value
)

trend_rows <- list()
for (cl in paste0("C", 1:4)) {
  successes <- as.numeric(tapply(df$cluster == cl, df$score_quartile, sum))
  totals <- as.numeric(tapply(rep(TRUE, nrow(df)), df$score_quartile, sum))
  pt <- prop.trend.test(successes, totals, score = seq_along(successes))
  trend_rows[[cl]] <- data.frame(
    target_cluster = cl,
    test = "Cochran-Armitage trend across score quartiles",
    statistic = unname(pt$statistic),
    p_value = pt$p.value,
    proportion_Q1 = successes[1] / totals[1],
    proportion_Q4 = successes[4] / totals[4],
    delta_Q4_minus_Q1 = successes[4] / totals[4] - successes[1] / totals[1],
    stringsAsFactors = FALSE
  )
}
trend_df <- rbindlist(trend_rows)

logistic_or <- rbindlist(lapply(paste0("C", 1:4), function(cl) or_from_glm(df, cl)))
logistic_or$BH_P <- p.adjust(logistic_or$P, method = "BH")

enrichment_top_bottom <- rbindlist(lapply(paste0("C", 1:4), function(cl) fisher_top_bottom(df, cl)))
enrichment_top_bottom$BH_P <- p.adjust(enrichment_top_bottom$P, method = "BH")
enrichment_top_bottom$Q1_proportion <- enrichment_top_bottom$Q1_target_n / enrichment_top_bottom$Q1_total_n
enrichment_top_bottom$Q4_proportion <- enrichment_top_bottom$Q4_target_n / enrichment_top_bottom$Q4_total_n
enrichment_top_bottom$Delta_Q4_minus_Q1 <- enrichment_top_bottom$Q4_proportion - enrichment_top_bottom$Q1_proportion

auc_rows <- list()
for (cl in paste0("C", 1:4)) {
  y <- as.integer(df$cluster == cl)
  auc <- auc_manual(y, df$Unified_ICB_score)
  auc_rows[[cl]] <- data.frame(
    target_cluster = cl,
    AUC_high_score_predicts_cluster = auc,
    OrientedAUC = max(auc, 1 - auc, na.rm = TRUE),
    Direction = ifelse(auc >= 0.5, "higher_score_more_cluster", "lower_score_more_cluster"),
    stringsAsFactors = FALSE
  )
}
auc_df <- rbindlist(auc_rows)

kw <- kruskal.test(Unified_ICB_score ~ cluster, data = df)
pairwise <- pairwise.wilcox.test(df$Unified_ICB_score, df$cluster, p.adjust.method = "BH")
pairwise_df <- as.data.frame(as.table(pairwise$p.value))
colnames(pairwise_df) <- c("Cluster1", "Cluster2", "BH_P")
pairwise_df <- pairwise_df[!is.na(pairwise_df$BH_P), ]

fwrite(df, file.path(out_dir, "TCGA_cluster4_class_samples_with_unified_ICB_score_quartile.txt"), sep = "\t")
fwrite(cluster_summary, file.path(out_dir, "cluster4_class_score_summary.txt"), sep = "\t")
fwrite(quartile_counts, file.path(out_dir, "score_quartile_cluster_composition.txt"), sep = "\t")
fwrite(chisq_df, file.path(out_dir, "cluster_composition_by_score_quartile_chisq.txt"), sep = "\t")
fwrite(trend_df, file.path(out_dir, "cluster_proportion_trend_across_score_quartiles.txt"), sep = "\t")
fwrite(logistic_or, file.path(out_dir, "cluster_membership_logistic_OR_per_1SD_score.txt"), sep = "\t")
fwrite(enrichment_top_bottom, file.path(out_dir, "cluster_enrichment_top_vs_bottom_score_quartile.txt"), sep = "\t")
fwrite(auc_df, file.path(out_dir, "score_predict_cluster_membership_auc.txt"), sep = "\t")
fwrite(data.frame(test = "Kruskal-Wallis score by cluster", statistic = unname(kw$statistic), df = unname(kw$parameter), p_value = kw$p.value),
       file.path(out_dir, "score_by_cluster_kruskal_test.txt"), sep = "\t")
fwrite(pairwise_df, file.path(out_dir, "score_by_cluster_pairwise_wilcox_BH.txt"), sep = "\t")

cluster_colors <- c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")

p1 <- ggplot(df, aes(x = cluster, y = Unified_ICB_score, fill = cluster)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black", linewidth = 0.28) +
  geom_jitter(width = 0.14, size = 0.35, alpha = 0.35) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 12) +
  labs(
    x = "cluster4_class",
    y = "Unified ICB Score",
    title = "Unified ICB Score by TCGA cluster4_class",
    subtitle = paste0("Kruskal-Wallis P = ", signif(kw$p.value, 3))
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "Unified_ICB_score_by_cluster4_class_boxplot.pdf"), p1, width = 5.2, height = 4.2)
ggsave(file.path(out_dir, "Unified_ICB_score_by_cluster4_class_boxplot.png"), p1, width = 5.2, height = 4.2, dpi = 300)

p2 <- ggplot(quartile_counts, aes(x = score_quartile, y = Proportion, fill = cluster)) +
  geom_col(color = "black", linewidth = 0.18, width = 0.72) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 11) +
  labs(
    x = "Unified ICB Score quartile",
    y = "Cluster proportion",
    fill = "cluster4_class",
    title = "Cluster composition shifts across Unified ICB Score quartiles",
    subtitle = paste0("Chi-square P = ", signif(chisq_cluster_quartile$p.value, 3))
  ) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cluster4_class_composition_by_score_quartile.pdf"), p2, width = 5.8, height = 4.2)
ggsave(file.path(out_dir, "cluster4_class_composition_by_score_quartile.png"), p2, width = 5.8, height = 4.2, dpi = 300)

p3 <- ggplot(trend_df, aes(x = target_cluster, y = delta_Q4_minus_Q1, fill = target_cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 11) +
  labs(
    x = "cluster4_class",
    y = "Proportion change: Q4 score - Q1 score",
    title = "Cluster enrichment from low to high Unified ICB Score"
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cluster4_class_Q4_minus_Q1_proportion_change.pdf"), p3, width = 5.2, height = 4)
ggsave(file.path(out_dir, "cluster4_class_Q4_minus_Q1_proportion_change.png"), p3, width = 5.2, height = 4, dpi = 300)

p4 <- ggplot(logistic_or, aes(x = target_cluster, y = OR_per_1SD_score, fill = target_cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.18, linewidth = 0.3) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dotted") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 11) +
  labs(
    x = "cluster4_class",
    y = "Odds ratio per 1 SD Unified ICB Score",
    title = "Association between Unified ICB Score and cluster membership"
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cluster4_class_logistic_OR_per_1SD_score.pdf"), p4, width = 5.2, height = 4)
ggsave(file.path(out_dir, "cluster4_class_logistic_OR_per_1SD_score.png"), p4, width = 5.2, height = 4, dpi = 300)

readme <- c(
  "Association between cluster4_class and Unified ICB Score",
  "Cluster file: <LOCAL_CLUSTER_ROOT>/cluster4_class.txt",
  "Score model: balanced ranger random forest using CLR-transformed 53 ImmuCellAI2 cell fractions.",
  "",
  "Cluster score summary:",
  capture.output(print(cluster_summary)),
  "",
  "Cluster composition by score quartile:",
  capture.output(print(quartile_counts)),
  "",
  "Trend of cluster proportions across score quartiles:",
  capture.output(print(trend_df)),
  "",
  "Logistic association per 1 SD score:",
  capture.output(print(logistic_or)),
  "",
  "Score predicting cluster membership AUC:",
  capture.output(print(auc_df))
)
writeLines(readme, file.path(out_dir, "README_cluster4_class_Unified_ICB_score_association.txt"))

print(cluster_summary)
print(quartile_counts)
print(chisq_df)
print(trend_df)
print(logistic_or)
print(enrichment_top_bottom)
print(auc_df)
