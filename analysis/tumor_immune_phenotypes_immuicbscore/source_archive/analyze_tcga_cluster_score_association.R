options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
out_dir <- file.path(fig4_dir, "TCGA_cluster_score_association")
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

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
})

set.seed(20260625)
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

message("Reading matrices...")
tcga_fraction <- read_fraction(tcga_fraction_file)
external_fraction <- read_fraction(external_fraction_file)
common_cells <- intersect(colnames(tcga_fraction), colnames(external_fraction))
tcga_fraction <- tcga_fraction[, common_cells, drop = FALSE]
external_fraction <- external_fraction[, common_cells, drop = FALSE]

message("Reading cluster labels...")
cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "cluster") %in% colnames(cluster_raw))) {
  stop("cluster4_class.txt must contain sample and cluster columns.")
}
cluster_df <- cluster_raw[, c("sample", "cluster")]
cluster_df$sample <- as.character(cluster_df$sample)
cluster_df$cluster <- paste0("Cluster", as.integer(cluster_df$cluster))
cluster_df <- cluster_df[cluster_df$sample %in% rownames(tcga_fraction), , drop = FALSE]
tcga_fraction <- tcga_fraction[cluster_df$sample, , drop = FALSE]

message("Training RF-CLR ICB response model on external ICB samples...")
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
  seed = 20260625,
  num.threads = n_threads
)
tcga_rf_icb_score <- as.numeric(
  predict(rf_fit, data = as.data.frame(scaled$test, check.names = FALSE))$predictions[, "1"]
)

message("Computing C3/C4 centroid-axis scores...")
tcga_clr_for_centroid <- clr_transform(tcga_fraction)
centroids_clr <- do.call(rbind, lapply(paste0("Cluster", 1:4), function(cl) {
  samples <- cluster_df$sample[cluster_df$cluster == cl]
  colMeans(tcga_clr_for_centroid[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_clr) <- paste0("Cluster", 1:4)

centroids_raw <- do.call(rbind, lapply(paste0("Cluster", 1:4), function(cl) {
  samples <- cluster_df$sample[cluster_df$cluster == cl]
  colMeans(tcga_fraction[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_raw) <- paste0("Cluster", 1:4)

c3_like_clr <- apply(tcga_clr_for_centroid, 1, safe_cor, y = centroids_clr["Cluster3", ], method = "pearson")
c4_like_clr <- apply(tcga_clr_for_centroid, 1, safe_cor, y = centroids_clr["Cluster4", ], method = "pearson")
c3_low_clr <- -c3_like_clr
c4_minus_c3_clr <- c4_like_clr - c3_like_clr
c4_minus_c3_raw <- apply(tcga_fraction, 1, safe_cor, y = centroids_raw["Cluster4", ], method = "spearman") -
  apply(tcga_fraction, 1, safe_cor, y = centroids_raw["Cluster3", ], method = "spearman")

score_df <- data.frame(
  sample = cluster_df$sample,
  cluster = factor(cluster_df$cluster, levels = paste0("Cluster", 1:4)),
  RF_CLR_ICB_response_score = tcga_rf_icb_score,
  C3_like_clr = c3_like_clr[cluster_df$sample],
  C3_low_clr = c3_low_clr[cluster_df$sample],
  C4_like_clr = c4_like_clr[cluster_df$sample],
  C4_minus_C3_clr = c4_minus_c3_clr[cluster_df$sample],
  C4_minus_C3_spearman_raw = c4_minus_c3_raw[cluster_df$sample],
  stringsAsFactors = FALSE
)

fwrite(score_df, file.path(out_dir, "TCGA_cluster_scores_sample_level.txt"), sep = "\t")

score_cols <- setdiff(colnames(score_df), c("sample", "cluster"))
summary_rows <- rbindlist(lapply(score_cols, function(sc) {
  dt <- as.data.table(score_df)
  dt[, .(
    N = .N,
    Mean = mean(get(sc), na.rm = TRUE),
    Median = median(get(sc), na.rm = TRUE),
    SD = sd(get(sc), na.rm = TRUE),
    Q25 = quantile(get(sc), 0.25, na.rm = TRUE),
    Q75 = quantile(get(sc), 0.75, na.rm = TRUE)
  ), by = cluster][, Score := sc][]
}), fill = TRUE)
setcolorder(summary_rows, c("Score", "cluster", "N", "Mean", "Median", "SD", "Q25", "Q75"))
fwrite(summary_rows, file.path(out_dir, "TCGA_cluster_score_summary.txt"), sep = "\t")

test_rows <- rbindlist(lapply(score_cols, function(sc) {
  kw <- kruskal.test(score_df[[sc]] ~ score_df$cluster)
  pair <- pairwise.wilcox.test(score_df[[sc]], score_df$cluster, p.adjust.method = "BH")
  get_pair_p <- function(a, b) {
    pv <- pair$p.value
    if (a %in% rownames(pv) && b %in% colnames(pv)) return(pv[a, b])
    if (b %in% rownames(pv) && a %in% colnames(pv)) return(pv[b, a])
    NA_real_
  }
  c4_vs_c3 <- get_pair_p("Cluster4", "Cluster3")
  c4_vs_c1 <- get_pair_p("Cluster4", "Cluster1")
  c4_vs_c2 <- get_pair_p("Cluster4", "Cluster2")
  data.frame(
    Score = sc,
    KruskalP = kw$p.value,
    BH_pairwise_P_C4_vs_C3 = c4_vs_c3,
    BH_pairwise_P_C4_vs_C1 = c4_vs_c1,
    BH_pairwise_P_C4_vs_C2 = c4_vs_c2,
    stringsAsFactors = FALSE
  )
}), fill = TRUE)
fwrite(test_rows, file.path(out_dir, "TCGA_cluster_score_stat_tests.txt"), sep = "\t")

cluster_order_check <- dcast(summary_rows, Score ~ cluster, value.var = "Mean")
cluster_order_check$HighestMeanCluster <- apply(
  cluster_order_check[, paste0("Cluster", 1:4), drop = FALSE],
  1,
  function(v) paste0("Cluster", which.max(v))
)
cluster_order_check$LowestMeanCluster <- apply(
  cluster_order_check[, paste0("Cluster", 1:4), drop = FALSE],
  1,
  function(v) paste0("Cluster", which.min(v))
)
fwrite(cluster_order_check, file.path(out_dir, "TCGA_cluster_score_highest_lowest_cluster.txt"), sep = "\t")

plot_df <- melt(
  as.data.table(score_df),
  id.vars = c("sample", "cluster"),
  measure.vars = score_cols,
  variable.name = "Score",
  value.name = "Value"
)

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F7E7E8", color = "black", linewidth = 0.25),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

p_all <- ggplot(plot_df, aes(x = cluster, y = Value, fill = cluster)) +
  geom_violin(width = 0.9, trim = TRUE, color = "grey25", linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, color = "black", linewidth = 0.25, fill = "white") +
  facet_wrap(~ Score, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF")) +
  labs(x = NULL, y = "Score", title = "Association between TCGA immune clusters and ICB-response scores") +
  theme_pub
ggsave(file.path(out_dir, "TCGA_cluster_score_association_all_scores.pdf"), p_all, width = 8.5, height = 7)
ggsave(file.path(out_dir, "TCGA_cluster_score_association_all_scores.png"), p_all, width = 8.5, height = 7, dpi = 300)

main_scores <- c("RF_CLR_ICB_response_score", "C3_low_clr", "C4_minus_C3_clr")
p_main <- ggplot(plot_df[Score %in% main_scores], aes(x = cluster, y = Value, fill = cluster)) +
  geom_violin(width = 0.9, trim = TRUE, color = "grey25", linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.3, fill = "white") +
  facet_wrap(~ Score, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF")) +
  labs(x = NULL, y = "Score", title = "TCGA cluster association with final ICB-response and C3/C4-axis scores") +
  theme_pub
ggsave(file.path(out_dir, "TCGA_cluster_score_association_main_scores.pdf"), p_main, width = 9, height = 3.6)
ggsave(file.path(out_dir, "TCGA_cluster_score_association_main_scores.png"), p_main, width = 9, height = 3.6, dpi = 300)

writeLines(c(
  "TCGA cluster-score association analysis",
  "RF_CLR_ICB_response_score: balanced ranger model trained on external ICB response labels using CLR-transformed 53 ImmuCellAI2 cell fractions, then projected onto TCGA samples.",
  "C3_low_clr: negative Pearson correlation between CLR-transformed sample fractions and TCGA Cluster3 centroid.",
  "C4_minus_C3_clr: C4-like CLR correlation minus C3-like CLR correlation.",
  "Use this analysis to show whether the cluster with better survival also carries higher predicted ICB-response/immune-active score."
), file.path(out_dir, "README_TCGA_cluster_score_association.txt"))

print(summary_rows)
print(test_rows)
print(cluster_order_check)
