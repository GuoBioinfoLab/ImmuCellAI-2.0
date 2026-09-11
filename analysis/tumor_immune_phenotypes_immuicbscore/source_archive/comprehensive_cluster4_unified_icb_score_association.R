options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_file <- "<LOCAL_CLUSTER_ROOT>/cluster4_class.txt"
score_file <- file.path(
  fig4_dir,
  "Unified_ICB_score_TCGA_cluster_association",
  "TCGA_samples_unified_ICB_score_by_cluster.txt"
)
survival_file <- file.path(fig4_dir, "survival_pan_all.RData")

out_dir <- file.path(fig4_dir, "cluster4_class_Unified_ICB_score_association_full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(grid)
})

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", hjust = 0),
      strip.background = element_rect(fill = "#F8E6E8", color = "black", linewidth = 0.35),
      strip.text = element_text(color = "black", face = "bold")
    )
}

cluster_cols <- c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")

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

jonckheere_stat <- function(x, group_ordered) {
  ok <- is.finite(x) & !is.na(group_ordered)
  x <- x[ok]
  group_ordered <- droplevels(group_ordered[ok])
  lev <- levels(group_ordered)
  j <- 0
  for (a in seq_len(length(lev) - 1L)) {
    xa <- x[group_ordered == lev[a]]
    for (b in (a + 1L):length(lev)) {
      xb <- sort(x[group_ordered == lev[b]])
      n_b <- length(xb)
      less_or_equal <- findInterval(xa, xb, rightmost.closed = TRUE)
      less_than <- findInterval(xa, xb, rightmost.closed = FALSE)
      equal_count <- less_or_equal - less_than
      greater_count <- n_b - less_or_equal
      j <- j + sum(greater_count + 0.5 * equal_count)
    }
  }
  as.numeric(j)
}

jonckheere_permutation <- function(x, group_ordered, n_perm = 5000L, seed = 20260626) {
  set.seed(seed)
  ok <- is.finite(x) & !is.na(group_ordered)
  x <- x[ok]
  group_ordered <- droplevels(group_ordered[ok])
  obs <- jonckheere_stat(x, group_ordered)
  perm <- numeric(n_perm)
  for (i in seq_len(n_perm)) perm[i] <- jonckheere_stat(x, sample(group_ordered))
  data.frame(
    JT_statistic = obs,
    permutation_mean = mean(perm),
    permutation_sd = sd(perm),
    z_permutation = (obs - mean(perm)) / sd(perm),
    p_increasing = (sum(perm >= obs) + 1) / (n_perm + 1),
    p_decreasing = (sum(perm <= obs) + 1) / (n_perm + 1),
    n_permutations = n_perm
  )
}

message("Reading cluster4_class and projected Unified ICB Score...")
cluster_raw <- fread(cluster_file, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "cluster") %in% colnames(cluster_raw))) {
  stop("cluster4_class.txt must contain columns named sample and cluster.")
}
cluster_df <- cluster_raw[, c("sample", "cluster")]
cluster_df$sample <- as.character(cluster_df$sample)
cluster_df$cluster <- paste0("C", as.integer(cluster_df$cluster))

score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
score_df <- score_df[, c("sample", "Unified_ICB_score")]
score_df$sample <- as.character(score_df$sample)

df <- merge(score_df, cluster_df, by = "sample")
df <- df[is.finite(df$Unified_ICB_score) & !is.na(df$cluster), ]
df$cluster <- factor(df$cluster, levels = paste0("C", 1:4))
df$score_z <- as.numeric(scale(df$Unified_ICB_score))

breaks <- quantile(df$Unified_ICB_score, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
breaks[1] <- breaks[1] - 1e-12
breaks[5] <- breaks[5] + 1e-12
df$score_quartile <- cut(
  df$Unified_ICB_score,
  breaks = breaks,
  labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"),
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

quartile_comp <- as.data.table(df)[, .N, by = .(score_quartile, cluster)]
quartile_comp[, Total := sum(N), by = score_quartile]
quartile_comp[, Proportion := N / Total]
quartile_comp <- quartile_comp[order(score_quartile, cluster)]

chisq_res <- suppressWarnings(chisq.test(table(df$score_quartile, df$cluster)))
chisq_df <- data.frame(
  test = "Cluster composition differs across Unified ICB Score quartiles",
  statistic = unname(chisq_res$statistic),
  df = unname(chisq_res$parameter),
  p_value = chisq_res$p.value
)

trend_rows <- list()
for (cl in paste0("C", 1:4)) {
  successes <- as.numeric(tapply(df$cluster == cl, df$score_quartile, sum))
  totals <- as.numeric(tapply(rep(TRUE, nrow(df)), df$score_quartile, sum))
  pt <- prop.trend.test(successes, totals, score = seq_along(successes))
  trend_rows[[cl]] <- data.frame(
    target_cluster = cl,
    statistic = unname(pt$statistic),
    p_value = pt$p.value,
    proportion_Q1 = successes[1] / totals[1],
    proportion_Q4 = successes[4] / totals[4],
    delta_Q4_minus_Q1 = successes[4] / totals[4] - successes[1] / totals[1],
    stringsAsFactors = FALSE
  )
}
trend_df <- rbindlist(trend_rows)
trend_df$BH_P <- p.adjust(trend_df$p_value, method = "BH")

logistic_or <- rbindlist(lapply(paste0("C", 1:4), function(cl) or_from_glm(df, cl)))
logistic_or$BH_P <- p.adjust(logistic_or$P, method = "BH")

top_bottom_rows <- list()
for (cl in paste0("C", 1:4)) {
  dd <- df[df$score_quartile %in% c("Q1 lowest", "Q4 highest"), ]
  tab <- table(HighScore = dd$score_quartile == "Q4 highest", TargetCluster = dd$cluster == cl)
  ft <- fisher.test(tab)
  top_bottom_rows[[cl]] <- data.frame(
    target_cluster = cl,
    OR_Q4_vs_Q1 = unname(ft$estimate),
    P = ft$p.value,
    Q1_proportion = mean(dd$cluster[dd$score_quartile == "Q1 lowest"] == cl),
    Q4_proportion = mean(dd$cluster[dd$score_quartile == "Q4 highest"] == cl),
    Delta_Q4_minus_Q1 = mean(dd$cluster[dd$score_quartile == "Q4 highest"] == cl) -
      mean(dd$cluster[dd$score_quartile == "Q1 lowest"] == cl),
    stringsAsFactors = FALSE
  )
}
top_bottom_df <- rbindlist(top_bottom_rows)
top_bottom_df$BH_P <- p.adjust(top_bottom_df$P, method = "BH")

auc_rows <- list()
for (cl in paste0("C", 1:4)) {
  auc <- auc_manual(as.integer(df$cluster == cl), df$Unified_ICB_score)
  auc_rows[[cl]] <- data.frame(
    target_cluster = cl,
    AUC_high_score_predicts_cluster = auc,
    Direction = ifelse(auc >= 0.5, "higher score enriched", "lower score enriched"),
    OrientedAUC = max(auc, 1 - auc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
auc_df <- rbindlist(auc_rows)

kw <- kruskal.test(Unified_ICB_score ~ cluster, data = df)
pairwise <- pairwise.wilcox.test(df$Unified_ICB_score, df$cluster, p.adjust.method = "BH")
pairwise_df <- as.data.frame(as.table(pairwise$p.value))
colnames(pairwise_df) <- c("Cluster1", "Cluster2", "BH_P")
pairwise_df <- pairwise_df[!is.na(pairwise_df$BH_P), ]

survival_outputs <- list()
if (file.exists(survival_file)) {
  load(survival_file)
  surv <- survival
  surv$sample <- as.character(surv$sample)
  surv_df <- merge(df, surv, by = "sample")
  surv_df <- surv_df[!is.na(surv_df$time) & !is.na(surv_df$status), ]

  surv_median <- do.call(rbind, lapply(levels(df$cluster), function(cl) {
    dd <- surv_df[surv_df$cluster == cl, ]
    sf <- survfit(Surv(time, status) ~ 1, data = dd)
    tab <- summary(sf)$table
    data.frame(
      cluster = cl,
      N_survival = nrow(dd),
      Events = sum(dd$status == 1),
      MedianSurvival = unname(tab["median"]),
      Lower95 = unname(tab["0.95LCL"]),
      Upper95 = unname(tab["0.95UCL"])
    )
  }))
  survival_order <- surv_median$cluster[order(surv_median$MedianSurvival)]
  df$survival_order <- factor(as.character(df$cluster), levels = survival_order, ordered = TRUE)
  jt <- jonckheere_permutation(df$Unified_ICB_score, df$survival_order, n_perm = 5000L)
  sp <- suppressWarnings(cor.test(df$Unified_ICB_score, as.integer(df$survival_order), method = "spearman", exact = FALSE))
  survival_trend <- cbind(
    cluster_order = paste(survival_order, collapse = " < "),
    jt,
    spearman_r = unname(sp$estimate),
    spearman_p = sp$p.value,
    stringsAsFactors = FALSE
  )
  survival_summary <- merge(cluster_summary, surv_median, by = "cluster", all.x = TRUE)
  survival_summary <- survival_summary[order(survival_summary$MedianSurvival), ]
  survival_outputs$survival_summary <- survival_summary
  survival_outputs$survival_trend <- survival_trend
}

fwrite(df, file.path(out_dir, "cluster4_class_samples_with_unified_ICB_score.txt"), sep = "\t")
fwrite(cluster_summary, file.path(out_dir, "cluster4_class_score_summary.txt"), sep = "\t")
fwrite(quartile_comp, file.path(out_dir, "cluster4_class_composition_by_score_quartile.txt"), sep = "\t")
fwrite(chisq_df, file.path(out_dir, "cluster4_class_score_quartile_chisq.txt"), sep = "\t")
fwrite(trend_df, file.path(out_dir, "cluster4_class_proportion_trend_across_score_quartiles.txt"), sep = "\t")
fwrite(logistic_or, file.path(out_dir, "cluster4_class_logistic_OR_per_1SD_score.txt"), sep = "\t")
fwrite(top_bottom_df, file.path(out_dir, "cluster4_class_top_vs_bottom_score_enrichment.txt"), sep = "\t")
fwrite(auc_df, file.path(out_dir, "unified_ICB_score_predict_cluster_membership_auc.txt"), sep = "\t")
fwrite(data.frame(test = "Kruskal-Wallis Unified ICB Score by cluster4_class", statistic = unname(kw$statistic), df = unname(kw$parameter), p_value = kw$p.value),
       file.path(out_dir, "cluster4_class_score_kruskal_test.txt"), sep = "\t")
fwrite(pairwise_df, file.path(out_dir, "cluster4_class_score_pairwise_wilcox_BH.txt"), sep = "\t")
if (length(survival_outputs)) {
  fwrite(survival_outputs$survival_summary, file.path(out_dir, "cluster4_class_score_and_survival_summary.txt"), sep = "\t")
  fwrite(survival_outputs$survival_trend, file.path(out_dir, "cluster4_class_score_survival_order_trend_test.txt"), sep = "\t")
}

p_box <- ggplot(df, aes(x = cluster, y = Unified_ICB_score, fill = cluster)) +
  geom_violin(width = 0.85, alpha = 0.55, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.36, outlier.shape = NA, color = "black", linewidth = 0.28, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.3, fill = "white", color = "black") +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "cluster4_class",
    y = "Unified ICB Score",
    title = "Unified ICB Score differs across TCGA cluster4_class",
    subtitle = paste0("Kruskal-Wallis P = ", signif(kw$p.value, 3))
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Figure_cluster4_score_violin_boxplot.pdf"), p_box, width = 5.2, height = 4.2)
ggsave(file.path(out_dir, "Figure_cluster4_score_violin_boxplot.png"), p_box, width = 5.2, height = 4.2, dpi = 300)

p_density <- ggplot(df, aes(x = Unified_ICB_score, color = cluster, fill = cluster)) +
  geom_density(alpha = 0.18, linewidth = 0.65) +
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(x = "Unified ICB Score", y = "Density", color = "cluster", fill = "cluster",
       title = "Score density by cluster4_class")
ggsave(file.path(out_dir, "Figure_cluster4_score_density.pdf"), p_density, width = 5.6, height = 4)
ggsave(file.path(out_dir, "Figure_cluster4_score_density.png"), p_density, width = 5.6, height = 4, dpi = 300)

p_stack <- ggplot(quartile_comp, aes(x = score_quartile, y = Proportion, fill = cluster)) +
  geom_col(color = "black", linewidth = 0.18, width = 0.72) +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "Unified ICB Score quartile",
    y = "Cluster proportion",
    fill = "cluster4_class",
    title = "cluster4_class composition shifts across score quartiles",
    subtitle = paste0("Chi-square P = ", signif(chisq_res$p.value, 3))
  )
ggsave(file.path(out_dir, "Figure_cluster4_composition_by_score_quartile.pdf"), p_stack, width = 5.8, height = 4.2)
ggsave(file.path(out_dir, "Figure_cluster4_composition_by_score_quartile.png"), p_stack, width = 5.8, height = 4.2, dpi = 300)

p_trend <- ggplot(quartile_comp, aes(x = score_quartile, y = Proportion, color = cluster, group = cluster)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(fill = cluster), shape = 21, color = "black", size = 2.4) +
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "Unified ICB Score quartile",
    y = "Cluster proportion",
    color = "cluster4_class",
    fill = "cluster4_class",
    title = "Cluster proportion trends along Unified ICB Score"
  )
ggsave(file.path(out_dir, "Figure_cluster4_proportion_trends_by_score_quartile.pdf"), p_trend, width = 5.8, height = 4.2)
ggsave(file.path(out_dir, "Figure_cluster4_proportion_trends_by_score_quartile.png"), p_trend, width = 5.8, height = 4.2, dpi = 300)

p_delta <- ggplot(trend_df, aes(x = target_cluster, y = delta_Q4_minus_Q1, fill = target_cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "cluster4_class",
    y = "Proportion change: Q4 - Q1",
    title = "High-score enrichment and depletion of cluster states"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Figure_cluster4_Q4_minus_Q1_enrichment.pdf"), p_delta, width = 5.2, height = 4)
ggsave(file.path(out_dir, "Figure_cluster4_Q4_minus_Q1_enrichment.png"), p_delta, width = 5.2, height = 4, dpi = 300)

p_or <- ggplot(logistic_or, aes(x = target_cluster, y = OR_per_1SD_score, fill = target_cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.18, linewidth = 0.3) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dotted") +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "cluster4_class",
    y = "Odds ratio per 1 SD score",
    title = "Unified ICB Score predicts cluster membership"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Figure_cluster4_logistic_OR_per_1SD_score.pdf"), p_or, width = 5.2, height = 4)
ggsave(file.path(out_dir, "Figure_cluster4_logistic_OR_per_1SD_score.png"), p_or, width = 5.2, height = 4, dpi = 300)

p_auc <- ggplot(auc_df, aes(x = target_cluster, y = OrientedAUC, fill = target_cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_hline(yintercept = 0.5, color = "grey40", linetype = "dotted") +
  scale_fill_manual(values = cluster_cols) +
  theme_pub() +
  labs(
    x = "cluster4_class",
    y = "Oriented AUC",
    title = "How well Unified ICB Score separates each cluster"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Figure_score_predict_cluster_membership_auc.pdf"), p_auc, width = 5.2, height = 4)
ggsave(file.path(out_dir, "Figure_score_predict_cluster_membership_auc.png"), p_auc, width = 5.2, height = 4, dpi = 300)

if (length(survival_outputs)) {
  surv_sum <- survival_outputs$survival_summary
  p_surv <- ggplot(surv_sum, aes(x = MedianSurvival, y = MeanScore, fill = cluster)) +
    geom_point(shape = 21, color = "black", size = 3.2) +
    geom_text(aes(label = cluster), vjust = -0.95, size = 3.4) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = cluster_cols) +
    theme_pub() +
    labs(
      x = "Median survival",
      y = "Mean Unified ICB Score",
      title = "Cluster-level score and survival association",
      subtitle = paste0(
        "Survival-order JT P = ",
        signif(survival_outputs$survival_trend$p_increasing, 3)
      )
    ) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "Figure_cluster_mean_score_vs_median_survival.pdf"), p_surv, width = 5.2, height = 4)
  ggsave(file.path(out_dir, "Figure_cluster_mean_score_vs_median_survival.png"), p_surv, width = 5.2, height = 4, dpi = 300)
}

pdf(file.path(out_dir, "Figure_cluster4_Unified_ICB_score_association_multipage.pdf"), width = 6.2, height = 4.8)
print(p_box)
print(p_density)
print(p_stack)
print(p_trend)
print(p_delta)
print(p_or)
print(p_auc)
if (exists("p_surv")) print(p_surv)
dev.off()

readme <- c(
  "Comprehensive association analysis between cluster4_class and Unified ICB Score",
  paste0("Cluster file: ", cluster_file),
  "Score model: balanced ranger random forest using CLR-transformed 53 ImmuCellAI2 cell fractions.",
  "",
  "Main interpretation:",
  "High Unified ICB Score is associated with depletion of C1 and enrichment of C4.",
  "",
  "Cluster summary:",
  capture.output(print(cluster_summary)),
  "",
  "Quartile composition:",
  capture.output(print(quartile_comp)),
  "",
  "Chi-square test:",
  capture.output(print(chisq_df)),
  "",
  "Cluster proportion trend across score quartiles:",
  capture.output(print(trend_df)),
  "",
  "Logistic OR per 1 SD score:",
  capture.output(print(logistic_or)),
  "",
  "Top-vs-bottom score enrichment:",
  capture.output(print(top_bottom_df)),
  "",
  "Score predicting cluster membership AUC:",
  capture.output(print(auc_df))
)
if (length(survival_outputs)) {
  readme <- c(
    readme,
    "",
    "Survival-linked cluster summary:",
    capture.output(print(survival_outputs$survival_summary)),
    "",
    "Survival-order trend test:",
    capture.output(print(survival_outputs$survival_trend))
  )
}
writeLines(readme, file.path(out_dir, "README_cluster4_class_Unified_ICB_score_association_full.txt"))

print(cluster_summary)
print(quartile_comp)
print(chisq_df)
print(trend_df)
print(logistic_or)
print(top_bottom_df)
print(auc_df)
if (length(survival_outputs)) print(survival_outputs$survival_trend)
