options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
out_dir <- file.path(fig4_dir, "TCGA_cluster_score_survival_association")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcga_fraction_file <- file.path(
  fig4_dir,
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000",
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
tcga_cluster_file <- file.path(cluster_dir, "cluster4_class.txt")
external_fraction_file <- file.path(
  base_dir,
  "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
external_meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
survival_file <- file.path(fig4_dir, "survival_pan_all.RData")

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
  library(survival)
  library(survminer)
})

set.seed(20260626)
n_threads <- min(8L, parallel::detectCores(logical = TRUE))

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

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train, na.rm = TRUE)
  scalev <- apply(x_train, 2, sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scalev, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scalev, "/")
  )
}

class_weights <- function(y) {
  ifelse(y == 1, 0.5 / mean(y == 1), 0.5 / mean(y == 0))
}

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

message("Reading ImmuCellAI2 matrices...")
tcga_fraction <- read_fraction(tcga_fraction_file)
external_fraction <- read_fraction(external_fraction_file)
common_cells <- intersect(colnames(tcga_fraction), colnames(external_fraction))
tcga_fraction <- tcga_fraction[, common_cells, drop = FALSE]
external_fraction <- external_fraction[, common_cells, drop = FALSE]

message("Reading cluster labels...")
cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
cluster_df <- cluster_raw[, c("sample", "cluster")]
cluster_df$sample <- as.character(cluster_df$sample)
cluster_df$cluster <- paste0("Cluster", as.integer(cluster_df$cluster))
cluster_df <- cluster_df[cluster_df$sample %in% rownames(tcga_fraction), , drop = FALSE]
tcga_fraction <- tcga_fraction[cluster_df$sample, , drop = FALSE]

message("Training best ICB ROC model: cell_clr53__ranger_balanced...")
external_meta <- fread(external_meta_file, data.table = FALSE, check.names = FALSE)
response <- rep(NA_integer_, nrow(external_meta))
response[external_meta$ResponseGroup == "Responder"] <- 1L
response[external_meta$ResponseGroup == "Non-responder"] <- 0L
response[is.na(response) & external_meta$Response == "R"] <- 1L
response[is.na(response) & external_meta$Response == "NR"] <- 0L
external_meta$response_binary <- response
external_meta <- external_meta[
  !is.na(external_meta$response_binary) & external_meta$Run %in% rownames(external_fraction),
  ,
  drop = FALSE
]
external_fraction <- external_fraction[external_meta$Run, , drop = FALSE]
y_icb <- external_meta$response_binary

external_clr <- clr_transform(external_fraction)
tcga_clr <- clr_transform(tcga_fraction)
colnames(external_clr) <- make.names(colnames(external_clr), unique = TRUE)
colnames(tcga_clr) <- make.names(colnames(tcga_clr), unique = TRUE)
scaled <- scale_train_test(external_clr, tcga_clr)

cw <- class_weights(y_icb)
rf_fit <- ranger(
  x = as.data.frame(scaled$train, check.names = FALSE),
  y = factor(y_icb, levels = c(0, 1)),
  probability = TRUE,
  classification = TRUE,
  num.trees = 1200,
  mtry = max(1, floor(sqrt(ncol(scaled$train)))),
  min.node.size = 8,
  class.weights = c("0" = cw[which(y_icb == 0)[1]], "1" = cw[which(y_icb == 1)[1]]),
  seed = 20260626,
  num.threads = n_threads
)
tcga_rf_icb_score <- as.numeric(
  predict(rf_fit, data = as.data.frame(scaled$test, check.names = FALSE))$predictions[, "1"]
)

message("Computing C3/C4 axis scores...")
tcga_clr_for_centroid <- clr_transform(tcga_fraction)
centroids_clr <- do.call(rbind, lapply(paste0("Cluster", 1:4), function(cl) {
  samples <- cluster_df$sample[cluster_df$cluster == cl]
  colMeans(tcga_clr_for_centroid[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_clr) <- paste0("Cluster", 1:4)

c3_like_clr <- apply(tcga_clr_for_centroid, 1, safe_cor, y = centroids_clr["Cluster3", ], method = "pearson")
c4_like_clr <- apply(tcga_clr_for_centroid, 1, safe_cor, y = centroids_clr["Cluster4", ], method = "pearson")
c3_low_clr <- -c3_like_clr
c4_minus_c3_clr <- c4_like_clr - c3_like_clr

score_df <- data.frame(
  sample = cluster_df$sample,
  patient = substr(cluster_df$sample, 1, 12),
  cluster = factor(cluster_df$cluster, levels = paste0("Cluster", 1:4)),
  RF_CLR_ICB_response_score = tcga_rf_icb_score,
  C3_like_clr = c3_like_clr[cluster_df$sample],
  C3_low_clr = c3_low_clr[cluster_df$sample],
  C4_like_clr = c4_like_clr[cluster_df$sample],
  C4_minus_C3_clr = c4_minus_c3_clr[cluster_df$sample],
  stringsAsFactors = FALSE
)

message("Reading survival...")
load(survival_file)
surv_df <- survival
surv_df$sample <- as.character(surv_df$sample)
surv_df$patient <- substr(surv_df$sample, 1, 12)
surv_df <- surv_df[order(surv_df$patient, surv_df$time, decreasing = TRUE), ]
surv_df <- surv_df[!duplicated(surv_df$patient), ]
merged <- merge(score_df, surv_df[, c("patient", "status", "time")], by = "patient", all.x = TRUE)
merged <- merged[is.finite(merged$time) & merged$time > 0 & !is.na(merged$status), , drop = FALSE]
merged$ScoreGroup <- ifelse(
  merged$RF_CLR_ICB_response_score >= median(merged$RF_CLR_ICB_response_score, na.rm = TRUE),
  "High RF-CLR score",
  "Low RF-CLR score"
)
merged$ScoreGroup <- factor(merged$ScoreGroup, levels = c("Low RF-CLR score", "High RF-CLR score"))

fwrite(merged, file.path(out_dir, "TCGA_cluster_RF_CLR_score_survival_sample_level.txt"), sep = "\t")

score_cols <- c("RF_CLR_ICB_response_score", "C3_low_clr", "C4_like_clr", "C4_minus_C3_clr")
summary_rows <- rbindlist(lapply(score_cols, function(sc) {
  dt <- as.data.table(merged)
  dt[, .(
    N = .N,
    Events = sum(status == 1, na.rm = TRUE),
    Mean = mean(get(sc), na.rm = TRUE),
    Median = median(get(sc), na.rm = TRUE),
    SD = sd(get(sc), na.rm = TRUE),
    Q25 = quantile(get(sc), 0.25, na.rm = TRUE),
    Q75 = quantile(get(sc), 0.75, na.rm = TRUE)
  ), by = cluster][, Score := sc][]
}), fill = TRUE)
setcolorder(summary_rows, c("Score", "cluster", "N", "Events", "Mean", "Median", "SD", "Q25", "Q75"))
fwrite(summary_rows, file.path(out_dir, "TCGA_cluster_score_summary_with_survival_samples.txt"), sep = "\t")

test_rows <- rbindlist(lapply(score_cols, function(sc) {
  kw <- kruskal.test(merged[[sc]] ~ merged$cluster)
  pair <- pairwise.wilcox.test(merged[[sc]], merged$cluster, p.adjust.method = "BH")
  get_pair_p <- function(a, b) {
    pv <- pair$p.value
    if (a %in% rownames(pv) && b %in% colnames(pv)) return(pv[a, b])
    if (b %in% rownames(pv) && a %in% colnames(pv)) return(pv[b, a])
    NA_real_
  }
  data.frame(
    Score = sc,
    KruskalP = kw$p.value,
    BH_P_Cluster4_vs_Cluster3 = get_pair_p("Cluster4", "Cluster3"),
    BH_P_Cluster4_vs_Cluster1 = get_pair_p("Cluster4", "Cluster1"),
    BH_P_Cluster4_vs_Cluster2 = get_pair_p("Cluster4", "Cluster2"),
    stringsAsFactors = FALSE
  )
}), fill = TRUE)
fwrite(test_rows, file.path(out_dir, "TCGA_cluster_score_stat_tests_with_survival_samples.txt"), sep = "\t")

highest_lowest <- dcast(summary_rows, Score ~ cluster, value.var = "Mean")
highest_lowest$HighestMeanCluster <- apply(highest_lowest[, paste0("Cluster", 1:4), drop = FALSE], 1, function(v) paste0("Cluster", which.max(v)))
highest_lowest$LowestMeanCluster <- apply(highest_lowest[, paste0("Cluster", 1:4), drop = FALSE], 1, function(v) paste0("Cluster", which.min(v)))
fwrite(highest_lowest, file.path(out_dir, "TCGA_cluster_score_highest_lowest_cluster.txt"), sep = "\t")

surv_by_cluster <- survfit(Surv(time, status) ~ cluster, data = merged)
surv_table <- as.data.frame(summary(surv_by_cluster)$table)
surv_tab <- data.frame(
  cluster = sub("^cluster=", "", rownames(surv_table)),
  N = as.numeric(surv_table[, "records"]),
  Events = as.numeric(surv_table[, "events"]),
  MedianSurvival = as.numeric(surv_table[, "median"]),
  Lower95 = as.numeric(surv_table[, "0.95LCL"]),
  Upper95 = as.numeric(surv_table[, "0.95UCL"]),
  row.names = NULL
)
fwrite(surv_tab, file.path(out_dir, "TCGA_cluster_survival_median_table.txt"), sep = "\t")

cox_rows <- rbind(
  tidy_cox(coxph(Surv(time, status) ~ relevel(cluster, ref = "Cluster4"), data = merged), "Cluster Cox, ref=Cluster4"),
  tidy_cox(coxph(Surv(time, status) ~ scale(RF_CLR_ICB_response_score), data = merged), "RF-CLR score continuous"),
  tidy_cox(coxph(Surv(time, status) ~ ScoreGroup, data = merged), "RF-CLR score high vs low"),
  tidy_cox(coxph(Surv(time, status) ~ scale(RF_CLR_ICB_response_score) + relevel(cluster, ref = "Cluster4"), data = merged), "RF-CLR score adjusted by cluster")
)
fwrite(cox_rows, file.path(out_dir, "TCGA_cluster_score_survival_cox_tests.txt"), sep = "\t")

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F7E7E8", color = "black", linewidth = 0.25),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

score_plot_df <- melt(
  as.data.table(merged),
  id.vars = c("sample", "cluster"),
  measure.vars = score_cols,
  variable.name = "Score",
  value.name = "Value"
)
p_scores <- ggplot(score_plot_df, aes(x = cluster, y = Value, fill = cluster)) +
  geom_violin(width = 0.9, trim = TRUE, color = "grey25", linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, color = "black", linewidth = 0.25, fill = "white") +
  facet_wrap(~ Score, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF")) +
  labs(x = NULL, y = "Score", title = "TCGA cluster association with ICB-response scores") +
  theme_pub
ggsave(file.path(out_dir, "TCGA_cluster_RF_CLR_and_C3C4_scores.pdf"), p_scores, width = 8.5, height = 6)
ggsave(file.path(out_dir, "TCGA_cluster_RF_CLR_and_C3C4_scores.png"), p_scores, width = 8.5, height = 6, dpi = 300)

p_main <- ggplot(merged, aes(x = cluster, y = RF_CLR_ICB_response_score, fill = cluster)) +
  geom_violin(width = 0.9, trim = TRUE, color = "grey25", linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, color = "black", linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF")) +
  labs(x = NULL, y = "Predicted ICB response score", title = "Best ROC score projected to TCGA immune clusters") +
  theme_pub
ggsave(file.path(out_dir, "TCGA_cluster_RF_CLR_ICB_response_score_main.pdf"), p_main, width = 5.2, height = 4)
ggsave(file.path(out_dir, "TCGA_cluster_RF_CLR_ICB_response_score_main.png"), p_main, width = 5.2, height = 4, dpi = 300)

p_km_cluster <- ggsurvplot(
  surv_by_cluster,
  data = merged,
  pval = TRUE,
  risk.table = FALSE,
  palette = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF"),
  legend.title = "Cluster",
  legend.labs = paste0("Cluster", 1:4)
)
ggsave(file.path(out_dir, "TCGA_cluster_KM_survival.pdf"), p_km_cluster$plot, width = 5.8, height = 4.6)
ggsave(file.path(out_dir, "TCGA_cluster_KM_survival.png"), p_km_cluster$plot, width = 5.8, height = 4.6, dpi = 300)

p_km_score <- ggsurvplot(
  survfit(Surv(time, status) ~ ScoreGroup, data = merged),
  data = merged,
  pval = TRUE,
  risk.table = FALSE,
  palette = c("#95A8AC", "#F66463"),
  legend.title = "RF-CLR score"
)
ggsave(file.path(out_dir, "TCGA_RF_CLR_score_high_low_KM_survival.pdf"), p_km_score$plot, width = 5.5, height = 4.4)
ggsave(file.path(out_dir, "TCGA_RF_CLR_score_high_low_KM_survival.png"), p_km_score$plot, width = 5.5, height = 4.4, dpi = 300)

writeLines(c(
  "TCGA cluster-score-survival association analysis",
  "RF_CLR_ICB_response_score was trained on external ICB response labels using CLR-transformed 53 ImmuCellAI2 cell fractions, then projected to TCGA samples.",
  "This analysis checks whether the survival-favorable TCGA cluster also has higher predicted ICB-response score.",
  paste0("Samples with survival used: ", nrow(merged))
), file.path(out_dir, "README_TCGA_cluster_score_survival_association.txt"))

print(summary_rows)
print(test_rows)
print(highest_lowest)
print(surv_tab)
print(cox_rows)
