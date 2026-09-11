options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
out_dir <- file.path(fig4_dir, "Unified_ICB_score_TCGA_cluster_association")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

icb_fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
icb_meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
tcga_fraction_file <- file.path(
  fig4_dir,
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000",
  "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
tcga_cluster_file <- file.path(cluster_dir, "cluster4_class.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
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

scale_train_test <- function(x_train, x_test) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  med <- apply(x_train, 2, median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(x_train))) {
    x_train[!is.finite(x_train[, j]), j] <- med[j]
    x_test[!is.finite(x_test[, j]), j] <- med[j]
  }
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 0
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
  center <- colMeans(x_train)
  scalev <- apply(x_train, 2, sd)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scalev, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scalev, "/")
  )
}

class_weights <- function(y) {
  ifelse(y == 1, 0.5 / mean(y == 1), 0.5 / mean(y == 0))
}

message("Reading ICB fractions and labels...")
icb_fraction <- read_fraction(icb_fraction_file)
icb_meta <- fread(icb_meta_file, data.table = FALSE, check.names = FALSE)

response <- rep(NA_integer_, nrow(icb_meta))
response[icb_meta$ResponseGroup == "Responder"] <- 1L
response[icb_meta$ResponseGroup == "Non-responder"] <- 0L
response[is.na(response) & icb_meta$Response == "R"] <- 1L
response[is.na(response) & icb_meta$Response == "NR"] <- 0L
icb_meta$response_binary <- response
icb_meta <- icb_meta[!is.na(icb_meta$response_binary) & icb_meta$Run %in% rownames(icb_fraction), , drop = FALSE]
icb_fraction <- icb_fraction[icb_meta$Run, , drop = FALSE]
y <- icb_meta$response_binary

message("Reading TCGA fractions and cluster labels...")
tcga_fraction <- read_fraction(tcga_fraction_file)
tcga_cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
tcga_cluster <- tcga_cluster_raw[, c("sample", "cluster")]
tcga_cluster$sample <- as.character(tcga_cluster$sample)
tcga_cluster$cluster <- paste0("C", as.integer(tcga_cluster$cluster))
tcga_cluster <- tcga_cluster[tcga_cluster$sample %in% rownames(tcga_fraction), , drop = FALSE]
tcga_fraction <- tcga_fraction[tcga_cluster$sample, , drop = FALSE]

common_cells <- intersect(colnames(icb_fraction), colnames(tcga_fraction))
icb_fraction <- icb_fraction[, common_cells, drop = FALSE]
tcga_fraction <- tcga_fraction[, common_cells, drop = FALSE]

message("Training unified ICB score model on all ICB samples...")
icb_clr <- clr_transform(icb_fraction)
tcga_clr <- clr_transform(tcga_fraction)
colnames(icb_clr) <- paste0("clr_", make.names(colnames(icb_fraction), unique = TRUE))
colnames(tcga_clr) <- colnames(icb_clr)

xs <- scale_train_test(icb_clr, tcga_clr)
cw <- class_weights(y)
fit <- ranger(
  x = as.data.frame(xs$train, check.names = FALSE),
  y = factor(y, levels = c(0, 1)),
  probability = TRUE,
  classification = TRUE,
  num.trees = 1200,
  mtry = max(1, floor(sqrt(ncol(xs$train)))),
  min.node.size = 8,
  class.weights = c("0" = cw[which(y == 0)[1]], "1" = cw[which(y == 1)[1]]),
  importance = "permutation",
  seed = 20260626,
  num.threads = n_threads
)

tcga_score <- as.numeric(predict(fit, data = as.data.frame(xs$test, check.names = FALSE))$predictions[, "1"])
icb_score_train <- as.numeric(predict(fit, data = as.data.frame(xs$train, check.names = FALSE))$predictions[, "1"])

tcga_score_df <- data.frame(
  sample = rownames(tcga_fraction),
  cluster = factor(tcga_cluster$cluster, levels = paste0("C", 1:4)),
  Unified_ICB_score = tcga_score,
  stringsAsFactors = FALSE
)

cluster_summary <- as.data.table(tcga_score_df)[, .(
  N = .N,
  MeanScore = mean(Unified_ICB_score, na.rm = TRUE),
  MedianScore = median(Unified_ICB_score, na.rm = TRUE),
  SDScore = sd(Unified_ICB_score, na.rm = TRUE),
  Q25 = quantile(Unified_ICB_score, 0.25, na.rm = TRUE),
  Q75 = quantile(Unified_ICB_score, 0.75, na.rm = TRUE)
), by = cluster][order(cluster)]

kw <- kruskal.test(Unified_ICB_score ~ cluster, data = tcga_score_df)
pairwise <- pairwise.wilcox.test(
  tcga_score_df$Unified_ICB_score,
  tcga_score_df$cluster,
  p.adjust.method = "BH"
)
pairwise_df <- as.data.frame(as.table(pairwise$p.value))
colnames(pairwise_df) <- c("Cluster1", "Cluster2", "BH_P")
pairwise_df <- pairwise_df[!is.na(pairwise_df$BH_P), ]

importance_df <- data.frame(
  feature = names(fit$variable.importance),
  importance = as.numeric(fit$variable.importance),
  row.names = NULL
)
importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]

icb_train_df <- data.frame(
  Run = icb_meta$Run,
  Response = icb_meta$Response,
  ResponseGroup = icb_meta$ResponseGroup,
  response_binary = y,
  Unified_ICB_score_train_all = icb_score_train,
  stringsAsFactors = FALSE
)

fwrite(tcga_score_df, file.path(out_dir, "TCGA_samples_unified_ICB_score_by_cluster.txt"), sep = "\t")
fwrite(cluster_summary, file.path(out_dir, "TCGA_cluster_unified_ICB_score_summary.txt"), sep = "\t")
fwrite(pairwise_df, file.path(out_dir, "TCGA_cluster_unified_ICB_score_pairwise_wilcox_BH.txt"), sep = "\t")
fwrite(data.frame(test = "Kruskal-Wallis", statistic = unname(kw$statistic), df = unname(kw$parameter), p_value = kw$p.value),
       file.path(out_dir, "TCGA_cluster_unified_ICB_score_kruskal_test.txt"), sep = "\t")
fwrite(importance_df, file.path(out_dir, "Unified_ICB_score_ranger_feature_importance.txt"), sep = "\t")
fwrite(icb_train_df, file.path(out_dir, "ICB_training_samples_unified_score_trained_all.txt"), sep = "\t")

p1 <- ggplot(tcga_score_df, aes(x = cluster, y = Unified_ICB_score, fill = cluster)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black", linewidth = 0.28) +
  geom_jitter(width = 0.14, size = 0.35, alpha = 0.35) +
  scale_fill_manual(values = c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")) +
  theme_bw(base_size = 12) +
  labs(
    x = "TCGA immune cluster",
    y = "Unified ICB score",
    title = "Projection of unified ICB response score onto TCGA immune clusters",
    subtitle = paste0("Kruskal-Wallis P = ", signif(kw$p.value, 3))
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "TCGA_cluster_unified_ICB_score_boxplot.pdf"), p1, width = 5.2, height = 4.2)
ggsave(file.path(out_dir, "TCGA_cluster_unified_ICB_score_boxplot.png"), p1, width = 5.2, height = 4.2, dpi = 300)

p2 <- ggplot(cluster_summary, aes(x = cluster, y = MeanScore, fill = cluster)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.65) +
  geom_errorbar(aes(ymin = MeanScore - SDScore, ymax = MeanScore + SDScore), width = 0.18, linewidth = 0.25) +
  scale_fill_manual(values = c("C1" = "#E15759FF", "C2" = "#4E79A7FF", "C3" = "#F28E2BFF", "C4" = "#76B7B2FF")) +
  theme_bw(base_size = 12) +
  labs(x = "TCGA immune cluster", y = "Mean Unified ICB score", title = "Mean score by TCGA cluster") +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "TCGA_cluster_unified_ICB_score_mean_barplot.pdf"), p2, width = 4.8, height = 3.8)
ggsave(file.path(out_dir, "TCGA_cluster_unified_ICB_score_mean_barplot.png"), p2, width = 4.8, height = 3.8, dpi = 300)

writeLines(c(
  "Unified ICB score projected onto TCGA C1-C4 immune clusters",
  "Score model: class-balanced ranger random forest trained on external ICB response labels.",
  "Predictors: CLR-transformed 53 ImmuCellAI2 cell fractions only.",
  "The trained score was projected back to TCGA immunotherapy-responsive tumor samples.",
  paste0("Kruskal-Wallis P: ", signif(kw$p.value, 6)),
  "",
  "Cluster score summary:",
  paste(capture.output(print(cluster_summary)), collapse = "\n")
), file.path(out_dir, "README_TCGA_cluster_unified_ICB_score_association.txt"))

print(cluster_summary)
print(data.frame(test = "Kruskal-Wallis", statistic = unname(kw$statistic), df = unname(kw$parameter), p_value = kw$p.value))
print(pairwise_df)
print(head(importance_df, 20))
