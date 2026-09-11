options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
score_dir <- file.path(fig4_dir, "Unified_ICB_score_TCGA_cluster_association")
survival_dir <- file.path(fig4_dir, "TCGA_cluster_score_survival_association")
out_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_trend_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(score_dir, "TCGA_samples_unified_ICB_score_by_cluster.txt")
surv_median_file <- file.path(survival_dir, "TCGA_cluster_survival_median_table.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

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
      # Count pairwise evidence that later/higher ordered groups have larger scores.
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
  for (i in seq_len(n_perm)) {
    perm[i] <- jonckheere_stat(x, sample(group_ordered))
  }
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

linear_trend_test <- function(df, order_col) {
  fit <- lm(Unified_ICB_score ~ .ord, data = data.frame(
    Unified_ICB_score = df$Unified_ICB_score,
    .ord = df[[order_col]]
  ))
  sm <- summary(fit)$coefficients
  data.frame(
    slope = sm[".ord", "Estimate"],
    t_value = sm[".ord", "t value"],
    p_value = sm[".ord", "Pr(>|t|)"],
    r_squared = summary(fit)$r.squared
  )
}

message("Reading score and survival order data...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
surv_df <- fread(surv_median_file, data.table = FALSE, check.names = FALSE)

surv_df$cluster_short <- sub("^Cluster", "C", surv_df$cluster)
survival_order <- surv_df$cluster_short[order(surv_df$MedianSurvival, decreasing = FALSE)]
survival_labels <- paste0(
  survival_order,
  " (median=", surv_df$MedianSurvival[match(survival_order, surv_df$cluster_short)], ")"
)

score_df$cluster <- as.character(score_df$cluster)
score_df$cluster_survival_order <- factor(score_df$cluster, levels = survival_order, ordered = TRUE)
score_df$survival_order_rank <- as.integer(score_df$cluster_survival_order)
score_df$cluster_numeric_order <- factor(score_df$cluster, levels = paste0("C", 1:4), ordered = TRUE)
score_df$numeric_order_rank <- as.integer(score_df$cluster_numeric_order)

mean_order <- as.data.table(score_df)[, .(MeanScore = mean(Unified_ICB_score)), by = cluster][order(MeanScore)]$cluster
score_df$cluster_score_mean_order <- factor(score_df$cluster, levels = mean_order, ordered = TRUE)
score_df$score_mean_order_rank <- as.integer(score_df$cluster_score_mean_order)

cluster_summary <- as.data.table(score_df)[, .(
  N = .N,
  MeanScore = mean(Unified_ICB_score, na.rm = TRUE),
  MedianScore = median(Unified_ICB_score, na.rm = TRUE),
  SDScore = sd(Unified_ICB_score, na.rm = TRUE),
  SurvivalMedian = surv_df$MedianSurvival[match(unique(cluster), surv_df$cluster_short)]
), by = cluster][order(SurvivalMedian)]

trend_rows <- list()

for (ord in list(
  survival_order = list(factor_col = "cluster_survival_order", rank_col = "survival_order_rank", order = survival_order),
  numeric_C1_to_C4 = list(factor_col = "cluster_numeric_order", rank_col = "numeric_order_rank", order = paste0("C", 1:4)),
  score_mean_order = list(factor_col = "cluster_score_mean_order", rank_col = "score_mean_order_rank", order = mean_order)
)) {
  name <- names(ord)
}

orders <- list(
  survival_order = list(factor_col = "cluster_survival_order", rank_col = "survival_order_rank", order = survival_order),
  numeric_C1_to_C4 = list(factor_col = "cluster_numeric_order", rank_col = "numeric_order_rank", order = paste0("C", 1:4)),
  score_mean_order = list(factor_col = "cluster_score_mean_order", rank_col = "score_mean_order_rank", order = mean_order)
)

trend_tables <- list()
for (nm in names(orders)) {
  ord <- orders[[nm]]
  jt <- jonckheere_permutation(score_df$Unified_ICB_score, score_df[[ord$factor_col]], n_perm = 5000L)
  lin <- linear_trend_test(score_df, ord$rank_col)
  sp <- suppressWarnings(cor.test(score_df$Unified_ICB_score, score_df[[ord$rank_col]], method = "spearman", exact = FALSE))
  kt <- suppressWarnings(cor.test(score_df$Unified_ICB_score, score_df[[ord$rank_col]], method = "kendall", exact = FALSE))
  trend_tables[[nm]] <- cbind(
    order_name = nm,
    cluster_order = paste(ord$order, collapse = " < "),
    jt,
    lin,
    spearman_r = unname(sp$estimate),
    spearman_p = sp$p.value,
    kendall_tau = unname(kt$estimate),
    kendall_p = kt$p.value,
    stringsAsFactors = FALSE
  )
}
trend_results <- rbindlist(trend_tables, fill = TRUE)

fwrite(score_df, file.path(out_dir, "TCGA_samples_unified_ICB_score_with_cluster_orders.txt"), sep = "\t")
fwrite(cluster_summary, file.path(out_dir, "TCGA_cluster_score_survival_order_summary.txt"), sep = "\t")
fwrite(trend_results, file.path(out_dir, "Unified_ICB_score_cluster_order_trend_tests.txt"), sep = "\t")

p1 <- ggplot(score_df, aes(x = cluster_survival_order, y = Unified_ICB_score, fill = cluster_survival_order)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black", linewidth = 0.28) +
  geom_jitter(width = 0.14, size = 0.35, alpha = 0.35) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "black") +
  scale_x_discrete(labels = survival_labels) +
  scale_fill_manual(values = c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")) +
  theme_bw(base_size = 11) +
  labs(
    x = "TCGA cluster ordered by median survival",
    y = "Unified ICB score",
    title = "Trend of Unified ICB score across survival-ordered TCGA clusters",
    subtitle = paste0(
      "Jonckheere increasing trend P = ",
      signif(trend_results$p_increasing[trend_results$order_name == "survival_order"], 3),
      "; Spearman r = ",
      signif(trend_results$spearman_r[trend_results$order_name == "survival_order"], 3)
    )
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "Unified_ICB_score_survival_ordered_cluster_trend_boxplot.pdf"), p1, width = 6.2, height = 4.5)
ggsave(file.path(out_dir, "Unified_ICB_score_survival_ordered_cluster_trend_boxplot.png"), p1, width = 6.2, height = 4.5, dpi = 300)

p2 <- ggplot(cluster_summary, aes(x = reorder(cluster, SurvivalMedian), y = MeanScore, group = 1)) +
  geom_line(color = "#379DA5", linewidth = 0.8) +
  geom_point(aes(fill = cluster), shape = 21, color = "black", size = 3) +
  geom_text(aes(label = paste0("median OS=", SurvivalMedian)), vjust = -0.9, size = 3) +
  scale_fill_manual(values = c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")) +
  theme_bw(base_size = 11) +
  labs(x = "Cluster ordered by median survival", y = "Mean Unified ICB score",
       title = "Mean Unified ICB score follows survival-favorable cluster order") +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "Unified_ICB_score_mean_trend_by_survival_order.pdf"), p2, width = 5.6, height = 4)
ggsave(file.path(out_dir, "Unified_ICB_score_mean_trend_by_survival_order.png"), p2, width = 5.6, height = 4, dpi = 300)

readme_lines <- c(
  "Unified ICB score trend test across TCGA clusters",
  paste0("Survival-order used as primary trend order: ", paste(survival_order, collapse = " < ")),
  paste0("Survival labels: ", paste(survival_labels, collapse = " < ")),
  "",
  "Cluster summary:",
  capture.output(print(cluster_summary)),
  "",
  "Trend tests:",
  capture.output(print(trend_results))
)
writeLines(readme_lines, file.path(out_dir, "README_Unified_ICB_score_cluster_trend_test.txt"))

print(cluster_summary)
print(trend_results)
