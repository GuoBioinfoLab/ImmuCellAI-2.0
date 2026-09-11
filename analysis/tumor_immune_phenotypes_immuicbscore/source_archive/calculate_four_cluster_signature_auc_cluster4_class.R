options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")

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
score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")

out_dir <- file.path(base_dir, "four_cluster_signature_response_auc_cluster4_class")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ranger)
  library(ggplot2)
})

set.seed(20260622)

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

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

cosine_sim <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!length(x) || den == 0) return(NA_real_)
  sum(x * y) / den
}

make_folds <- function(y, k = 5) {
  folds <- integer(length(y))
  for (yy in sort(unique(y))) {
    idx <- sample(which(y == yy))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
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

cv_glmnet_auc <- function(x, y, alpha = 0, n_repeats = 30, k = 5) {
  rows <- vector("list", n_repeats)
  preds <- vector("list", n_repeats)
  for (rr in seq_len(n_repeats)) {
    folds <- make_folds(y, k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      train <- which(folds != ff)
      test <- which(folds == ff)
      xs <- scale_train_test(as.matrix(x[train, , drop = FALSE]), as.matrix(x[test, , drop = FALSE]))
      fit <- suppressWarnings(cv.glmnet(
        xs$train,
        y[train],
        family = "binomial",
        alpha = alpha,
        type.measure = "auc",
        nfolds = 5
      ))
      pred[test] <- as.numeric(predict(fit, newx = xs$test, s = "lambda.min", type = "response"))
    }
    rows[[rr]] <- data.frame(
      repeat_id = rr,
      AUC = auc_manual(y, pred),
      Mean_R = mean(pred[y == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y == 0], na.rm = TRUE)
    )
    preds[[rr]] <- data.frame(repeat_id = rr, sample = names(y), y = y, prediction = pred)
  }
  list(auc = do.call(rbind, rows), pred = do.call(rbind, preds))
}

cv_ranger_auc <- function(x, y, n_repeats = 30, k = 5) {
  rows <- vector("list", n_repeats)
  preds <- vector("list", n_repeats)
  for (rr in seq_len(n_repeats)) {
    folds <- make_folds(y, k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      train <- which(folds != ff)
      test <- which(folds == ff)
      xs <- scale_train_test(as.matrix(x[train, , drop = FALSE]), as.matrix(x[test, , drop = FALSE]))
      dat <- data.frame(y = factor(y[train], levels = c(0, 1)), as.data.frame(xs$train, check.names = FALSE))
      fit <- ranger(
        y ~ .,
        data = dat,
        probability = TRUE,
        classification = TRUE,
        num.trees = 800,
        mtry = max(1, floor(sqrt(ncol(xs$train)))),
        min.node.size = 10,
        seed = sample.int(1e6, 1),
        num.threads = min(8L, parallel::detectCores(logical = TRUE))
      )
      pred[test] <- predict(fit, data = as.data.frame(xs$test, check.names = FALSE))$predictions[, "1"]
    }
    rows[[rr]] <- data.frame(
      repeat_id = rr,
      AUC = auc_manual(y, pred),
      Mean_R = mean(pred[y == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y == 0], na.rm = TRUE)
    )
    preds[[rr]] <- data.frame(repeat_id = rr, sample = names(y), y = y, prediction = pred)
  }
  list(auc = do.call(rbind, rows), pred = do.call(rbind, preds))
}

read_fraction <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  rownames(x) <- as.character(x[[1]])
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

message("Reading TCGA and external ImmuCellAI2 fractions...")
tcga_fraction <- read_fraction(tcga_fraction_file)
external_fraction <- read_fraction(external_fraction_file)
common_cells <- intersect(colnames(tcga_fraction), colnames(external_fraction))
tcga_fraction <- tcga_fraction[, common_cells, drop = FALSE]
external_fraction <- external_fraction[, common_cells, drop = FALSE]

message("Reading cluster4_class.txt...")
tcga_cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "cluster") %in% colnames(tcga_cluster_raw))) {
  stop("cluster4_class.txt must contain columns named sample and cluster.")
}
tcga_cluster <- tcga_cluster_raw[, c("sample", "cluster")]
tcga_cluster$sample <- as.character(tcga_cluster$sample)
tcga_cluster$cluster <- paste0("C", as.integer(tcga_cluster$cluster))
tcga_cluster <- tcga_cluster[tcga_cluster$sample %in% rownames(tcga_fraction), , drop = FALSE]

cluster_sizes <- as.data.frame(table(tcga_cluster$cluster), stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("cluster", "n_samples")
fwrite(cluster_sizes, file.path(out_dir, "cluster4_class_sizes_used_for_centroids.txt"), sep = "\t")

centroids <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- tcga_cluster$sample[tcga_cluster$cluster == cl]
  colMeans(tcga_fraction[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids) <- paste0("C", 1:4)
fwrite(data.frame(cluster = rownames(centroids), centroids, check.names = FALSE),
       file.path(out_dir, "TCGA_cluster4_class_centroids_53celltypes.txt"), sep = "\t")

score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
response_binary <- rep(NA_integer_, nrow(score_df))
response_binary[score_df$ResponseGroup == "Responder"] <- 1L
response_binary[score_df$ResponseGroup == "Non-responder"] <- 0L
response_binary[is.na(response_binary) & score_df$Response == "R"] <- 1L
response_binary[is.na(response_binary) & score_df$Response == "NR"] <- 0L
score_df$response_binary <- response_binary
score_df <- score_df[!is.na(score_df$response_binary) & score_df$Run %in% rownames(external_fraction), , drop = FALSE]
external_fraction <- external_fraction[score_df$Run, , drop = FALSE]
y <- score_df$response_binary
names(y) <- score_df$Run

message("External matched samples: ", length(y), "; responders: ", sum(y == 1), "; nonresponders: ", sum(y == 0))

cluster_features <- data.frame(sample = rownames(external_fraction), stringsAsFactors = FALSE)
for (cl in rownames(centroids)) {
  cluster_features[[paste0(cl, "_spearman")]] <- apply(external_fraction, 1, safe_cor, y = centroids[cl, ], method = "spearman")
  cluster_features[[paste0(cl, "_pearson")]] <- apply(external_fraction, 1, safe_cor, y = centroids[cl, ], method = "pearson")
  cluster_features[[paste0(cl, "_cosine")]] <- apply(external_fraction, 1, cosine_sim, y = centroids[cl, ])
}
for (a in paste0("C", 1:4)) {
  for (b in paste0("C", 1:4)) {
    if (a < b) {
      cluster_features[[paste0(b, "_minus_", a, "_spearman")]] <-
        cluster_features[[paste0(b, "_spearman")]] - cluster_features[[paste0(a, "_spearman")]]
    }
  }
}
fwrite(cluster_features, file.path(out_dir, "external_four_cluster_similarity_features_cluster4_class.txt"), sep = "\t")

univariate <- data.frame(feature = setdiff(colnames(cluster_features), "sample"), AUC = NA_real_)
for (i in seq_len(nrow(univariate))) {
  a <- auc_manual(y, cluster_features[[univariate$feature[i]]])
  univariate$AUC[i] <- max(a, 1 - a, na.rm = TRUE)
}
univariate <- univariate[order(univariate$AUC, decreasing = TRUE), ]

feature_sets <- list(
  four_spearman = paste0("C", 1:4, "_spearman"),
  four_pearson = paste0("C", 1:4, "_pearson"),
  four_cosine = paste0("C", 1:4, "_cosine"),
  all12_similarity = grep("_(spearman|pearson|cosine)$", colnames(cluster_features), value = TRUE),
  four_spearman_plus_deltas = c(paste0("C", 1:4, "_spearman"), grep("_minus_.*_spearman$", colnames(cluster_features), value = TRUE))
)

cv_summary <- list()
cv_repeat <- list()
cv_predictions <- list()
for (nm in names(feature_sets)) {
  x <- cluster_features[, feature_sets[[nm]], drop = FALSE]
  for (method in c("glmnet_ridge", "ranger")) {
    message("Running ", nm, " / ", method)
    res <- if (method == "glmnet_ridge") {
      cv_glmnet_auc(x, y, alpha = 0, n_repeats = 30, k = 5)
    } else {
      cv_ranger_auc(x, y, n_repeats = 30, k = 5)
    }
    model_name <- paste(nm, method, sep = "__")
    cv_repeat[[model_name]] <- transform(res$auc, model = model_name)
    cv_predictions[[model_name]] <- transform(res$pred, model = model_name)
    cv_summary[[model_name]] <- data.frame(
      model = model_name,
      feature_set = nm,
      method = method,
      MeanAUC = mean(res$auc$AUC, na.rm = TRUE),
      MedianAUC = median(res$auc$AUC, na.rm = TRUE),
      SDAUC = sd(res$auc$AUC, na.rm = TRUE),
      MaxAUC = max(res$auc$AUC, na.rm = TRUE),
      MinAUC = min(res$auc$AUC, na.rm = TRUE),
      Mean_R = mean(res$auc$Mean_R, na.rm = TRUE),
      Mean_NR = mean(res$auc$Mean_NR, na.rm = TRUE)
    )
  }
}

cv_repeat <- do.call(rbind, cv_repeat)
cv_predictions <- do.call(rbind, cv_predictions)
cv_summary <- do.call(rbind, cv_summary)
cv_summary <- cv_summary[order(cv_summary$MeanAUC, decreasing = TRUE), ]

fwrite(univariate, file.path(out_dir, "four_cluster_univariate_auc_cluster4_class.txt"), sep = "\t")
fwrite(cv_repeat, file.path(out_dir, "four_cluster_repeated5fold_auc_by_repeat_cluster4_class.txt"), sep = "\t")
fwrite(cv_predictions, file.path(out_dir, "four_cluster_repeated5fold_predictions_cluster4_class.txt"), sep = "\t")
fwrite(cv_summary, file.path(out_dir, "four_cluster_repeated5fold_auc_summary_cluster4_class.txt"), sep = "\t")

p <- ggplot(cv_summary, aes(x = reorder(model, MeanAUC), y = MeanAUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Mean repeated 5-fold CV AUC",
       title = "Four-cluster features from cluster4_class.txt")
ggsave(file.path(out_dir, "four_cluster_feature_auc_barplot_cluster4_class.pdf"), p, width = 7.2, height = 4.8)
ggsave(file.path(out_dir, "four_cluster_feature_auc_barplot_cluster4_class.png"), p, width = 7.2, height = 4.8, dpi = 300)

print(cluster_sizes)
print(cv_summary)
print(head(univariate, 20))
