options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
score_dir <- file.path(fig4_dir, "Unified_ICB_score_TCGA_cluster_association")
out_dir <- file.path(fig4_dir, "Unified_ICB_score_original_cluster_trend_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(score_dir, "TCGA_samples_unified_ICB_score_by_cluster.txt")
original_cluster_file <- file.path(cluster_dir, "cluster_TCGA.k=4.consensusClass.csv")
survival_file <- file.path(fig4_dir, "survival_pan_all.RData")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
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

linear_trend_test <- function(df, rank_col) {
  fit <- lm(Unified_ICB_score ~ .ord, data = data.frame(
    Unified_ICB_score = df$Unified_ICB_score,
    .ord = df[[rank_col]]
  ))
  sm <- summary(fit)$coefficients
  data.frame(
    slope = sm[".ord", "Estimate"],
    t_value = sm[".ord", "t value"],
    p_value = sm[".ord", "Pr(>|t|)"],
    r_squared = summary(fit)$r.squared
  )
}

message("Reading score, original cluster, and survival data...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
score_df <- score_df[, c("sample", "Unified_ICB_score")]

orig_cluster <- fread(original_cluster_file, header = FALSE, data.table = FALSE)
colnames(orig_cluster) <- c("sample", "cluster_id")
orig_cluster$sample <- gsub('"', "", as.character(orig_cluster$sample))
orig_cluster$cluster <- paste0("Cluster", as.integer(orig_cluster$cluster_id))

load(survival_file)
surv <- survival
surv$sample <- as.character(surv$sample)

df <- merge(score_df, orig_cluster[, c("sample", "cluster")], by = "sample")
df <- merge(df, surv, by = "sample", all.x = TRUE)
df$cluster <- factor(df$cluster, levels = paste0("Cluster", 1:4))

surv_df <- df[!is.na(df$time) & !is.na(df$status), ]
surv_median <- do.call(rbind, lapply(levels(df$cluster), function(cl) {
  dd <- surv_df[surv_df$cluster == cl, ]
  sf <- survfit(Surv(time, status) ~ 1, data = dd)
  med <- summary(sf)$table
  data.frame(
    cluster = cl,
    N_survival = nrow(dd),
    Events = sum(dd$status == 1, na.rm = TRUE),
    MedianSurvival = unname(med["median"]),
    Lower95 = unname(med["0.95LCL"]),
    Upper95 = unname(med["0.95UCL"])
  )
}))

survival_order <- surv_median$cluster[order(surv_median$MedianSurvival, decreasing = FALSE)]
score_mean_order <- as.data.table(df)[, .(MeanScore = mean(Unified_ICB_score)), by = cluster][order(MeanScore)]$cluster

df$cluster_survival_order <- factor(df$cluster, levels = survival_order, ordered = TRUE)
df$survival_order_rank <- as.integer(df$cluster_survival_order)
df$cluster_numeric_order <- factor(df$cluster, levels = paste0("Cluster", 1:4), ordered = TRUE)
df$numeric_order_rank <- as.integer(df$cluster_numeric_order)
df$cluster_score_mean_order <- factor(df$cluster, levels = score_mean_order, ordered = TRUE)
df$score_mean_order_rank <- as.integer(df$cluster_score_mean_order)

cluster_summary <- as.data.table(df)[, .(
  N = .N,
  MeanScore = mean(Unified_ICB_score, na.rm = TRUE),
  MedianScore = median(Unified_ICB_score, na.rm = TRUE),
  SDScore = sd(Unified_ICB_score, na.rm = TRUE),
  N_survival = sum(!is.na(time) & !is.na(status))
), by = cluster]
cluster_summary <- merge(cluster_summary, surv_median, by = "cluster", all.x = TRUE)
cluster_summary <- cluster_summary[order(MedianSurvival), ]

orders <- list(
  survival_order = list(factor_col = "cluster_survival_order", rank_col = "survival_order_rank", order = survival_order),
  numeric_Cluster1_to_Cluster4 = list(factor_col = "cluster_numeric_order", rank_col = "numeric_order_rank", order = paste0("Cluster", 1:4)),
  score_mean_order = list(factor_col = "cluster_score_mean_order", rank_col = "score_mean_order_rank", order = score_mean_order)
)

trend_tables <- list()
for (nm in names(orders)) {
  ord <- orders[[nm]]
  jt <- jonckheere_permutation(df$Unified_ICB_score, df[[ord$factor_col]], n_perm = 5000L)
  lin <- linear_trend_test(df, ord$rank_col)
  sp <- suppressWarnings(cor.test(df$Unified_ICB_score, df[[ord$rank_col]], method = "spearman", exact = FALSE))
  kt <- suppressWarnings(cor.test(df$Unified_ICB_score, df[[ord$rank_col]], method = "kendall", exact = FALSE))
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

cox_cluster <- coxph(Surv(time, status) ~ relevel(cluster, ref = survival_order[length(survival_order)]), data = surv_df)
cox_score <- coxph(Surv(time, status) ~ scale(Unified_ICB_score), data = surv_df)
cox_score_cluster <- coxph(Surv(time, status) ~ scale(Unified_ICB_score) + relevel(cluster, ref = survival_order[length(survival_order)]), data = surv_df)

tidy_cox <- function(fit, model_name) {
  sm <- summary(fit)
  data.frame(
    Model = model_name,
    Term = rownames(sm$coefficients),
    HR = sm$coefficients[, "exp(coef)"],
    Lower95 = sm$conf.int[, "lower .95"],
    Upper95 = sm$conf.int[, "upper .95"],
    P = sm$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}
cox_tests <- rbind(
  tidy_cox(cox_cluster, paste0("Original cluster Cox, ref=", survival_order[length(survival_order)])),
  tidy_cox(cox_score, "Unified ICB score continuous"),
  tidy_cox(cox_score_cluster, paste0("Unified ICB score adjusted by original cluster, ref=", survival_order[length(survival_order)]))
)

fwrite(df, file.path(out_dir, "TCGA_samples_unified_ICB_score_original_cluster_orders.txt"), sep = "\t")
fwrite(cluster_summary, file.path(out_dir, "original_cluster_score_survival_order_summary.txt"), sep = "\t")
fwrite(trend_results, file.path(out_dir, "original_cluster_unified_ICB_score_trend_tests.txt"), sep = "\t")
fwrite(cox_tests, file.path(out_dir, "original_cluster_unified_ICB_score_cox_tests.txt"), sep = "\t")

survival_labels <- paste0(
  survival_order,
  " (median=", surv_median$MedianSurvival[match(survival_order, surv_median$cluster)], ")"
)

p1 <- ggplot(df, aes(x = cluster_survival_order, y = Unified_ICB_score, fill = cluster_survival_order)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black", linewidth = 0.28) +
  geom_jitter(width = 0.14, size = 0.35, alpha = 0.35) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "black") +
  scale_x_discrete(labels = survival_labels) +
  scale_fill_manual(values = c("Cluster1" = "#E15759FF", "Cluster2" = "#4E79A7FF", "Cluster3" = "#F28E2BFF", "Cluster4" = "#76B7B2FF")) +
  theme_bw(base_size = 11) +
  labs(
    x = "Original TCGA cluster ordered by median survival",
    y = "Unified ICB score",
    title = "Unified ICB score trend across original survival-ordered clusters",
    subtitle = paste0(
      "JT increasing P = ",
      signif(trend_results$p_increasing[trend_results$order_name == "survival_order"], 3),
      "; Spearman r = ",
      signif(trend_results$spearman_r[trend_results$order_name == "survival_order"], 3)
    )
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "original_cluster_unified_ICB_score_survival_order_trend_boxplot.pdf"), p1, width = 6.4, height = 4.5)
ggsave(file.path(out_dir, "original_cluster_unified_ICB_score_survival_order_trend_boxplot.png"), p1, width = 6.4, height = 4.5, dpi = 300)

p2 <- ggplot(cluster_summary, aes(x = reorder(cluster, MedianSurvival), y = MeanScore, group = 1)) +
  geom_line(color = "#379DA5", linewidth = 0.8) +
  geom_point(aes(fill = cluster), shape = 21, color = "black", size = 3) +
  geom_text(aes(label = paste0("median=", MedianSurvival)), vjust = -0.9, size = 3) +
  scale_fill_manual(values = c("Cluster1" = "#E15759FF", "Cluster2" = "#4E79A7FF", "Cluster3" = "#F28E2BFF", "Cluster4" = "#76B7B2FF")) +
  theme_bw(base_size = 11) +
  labs(x = "Original cluster ordered by median survival", y = "Mean Unified ICB score",
       title = "Mean Unified ICB score across original survival-ordered clusters") +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "original_cluster_unified_ICB_score_mean_trend_by_survival_order.pdf"), p2, width = 5.8, height = 4)
ggsave(file.path(out_dir, "original_cluster_unified_ICB_score_mean_trend_by_survival_order.png"), p2, width = 5.8, height = 4, dpi = 300)

writeLines(c(
  "Unified ICB score trend test across original ConsensusClusterPlus clusters",
  paste0("Original cluster file: ", original_cluster_file),
  paste0("Primary survival order: ", paste(survival_order, collapse = " < ")),
  "",
  "Cluster summary:",
  capture.output(print(cluster_summary)),
  "",
  "Trend tests:",
  capture.output(print(trend_results)),
  "",
  "Cox tests:",
  capture.output(print(cox_tests))
), file.path(out_dir, "README_original_cluster_unified_ICB_score_trend_test.txt"))

print(cluster_summary)
print(trend_results)
print(cox_tests)
