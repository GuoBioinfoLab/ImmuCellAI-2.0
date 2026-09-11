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

out_dir <- file.path(base_dir, "C3_C4_axis_ICB_response_scores")
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

score_summary <- function(y, score) {
  a <- auc_manual(y, score)
  data.frame(
    AUC = a,
    OrientedAUC = max(a, 1 - a, na.rm = TRUE),
    Direction = ifelse(is.na(a), NA_character_, ifelse(a >= 0.5, "higher_in_R", "higher_in_NR")),
    Mean_R = mean(score[y == 1], na.rm = TRUE),
    Mean_NR = mean(score[y == 0], na.rm = TRUE),
    Median_R = median(score[y == 1], na.rm = TRUE),
    Median_NR = median(score[y == 0], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
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

row_mean_existing <- function(mat, cells) {
  cells <- intersect(cells, colnames(mat))
  if (!length(cells)) return(rep(NA_real_, nrow(mat)))
  rowMeans(mat[, cells, drop = FALSE], na.rm = TRUE)
}

row_sum_existing <- function(mat, cells) {
  cells <- intersect(cells, colnames(mat))
  if (!length(cells)) return(rep(NA_real_, nrow(mat)))
  rowSums(mat[, cells, drop = FALSE], na.rm = TRUE)
}

read_fraction <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  rownames(x) <- as.character(x[[1]])
  x[[1]] <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  m
}

asin_transform <- function(x) asin(sqrt(pmax(x, 0)))

clr_transform <- function(x, eps = 1e-5) {
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx), "-")
}

z_by_train <- function(train, test) {
  center <- colMeans(train, na.rm = TRUE)
  scalev <- apply(train, 2, sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(
    train = sweep(sweep(train, 2, center, "-"), 2, scalev, "/"),
    test = sweep(sweep(test, 2, center, "-"), 2, scalev, "/")
  )
}

message("Reading matrices...")
tcga_fraction <- read_fraction(tcga_fraction_file)
external_fraction <- read_fraction(external_fraction_file)
common_cells <- intersect(colnames(tcga_fraction), colnames(external_fraction))
tcga_fraction <- tcga_fraction[, common_cells, drop = FALSE]
external_fraction <- external_fraction[, common_cells, drop = FALSE]

cluster_raw <- fread(tcga_cluster_file, data.table = FALSE, check.names = FALSE)
cluster <- cluster_raw[, c("sample", "cluster")]
cluster$sample <- as.character(cluster$sample)
cluster$cluster <- paste0("C", as.integer(cluster$cluster))
cluster <- cluster[cluster$sample %in% rownames(tcga_fraction), , drop = FALSE]

cluster_sizes <- as.data.frame(table(cluster$cluster), stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("cluster", "n_samples")
fwrite(cluster_sizes, file.path(out_dir, "cluster4_class_sizes_used.txt"), sep = "\t")

centroids <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- cluster$sample[cluster$cluster == cl]
  colMeans(tcga_fraction[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids) <- paste0("C", 1:4)
fwrite(data.frame(cluster = rownames(centroids), centroids, check.names = FALSE),
       file.path(out_dir, "cluster4_class_centroids.txt"), sep = "\t")

c4_minus_c3 <- centroids["C4", ] - centroids["C3", ]
c3_minus_c4 <- -c4_minus_c3
c4_high <- names(sort(c4_minus_c3, decreasing = TRUE))
c3_high <- names(sort(c3_minus_c4, decreasing = TRUE))

top_table <- data.frame(
  Rank = seq_len(length(c4_high)),
  C4_high_cell = c4_high,
  C4_minus_C3 = c4_minus_c3[c4_high],
  C3_high_cell = c3_high,
  C3_minus_C4 = c3_minus_c4[c3_high],
  row.names = NULL,
  check.names = FALSE
)
fwrite(top_table, file.path(out_dir, "C4_high_and_C3_high_cells_ranked.txt"), sep = "\t")

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
message("External matched samples: ", length(y), "; R=", sum(y == 1), "; NR=", sum(y == 0))

external_asin <- asin_transform(external_fraction)
external_clr <- clr_transform(external_fraction)
tcga_asin <- asin_transform(tcga_fraction)
tcga_clr <- clr_transform(tcga_fraction)

centroids_asin <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- cluster$sample[cluster$cluster == cl]
  colMeans(tcga_asin[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_asin) <- paste0("C", 1:4)

centroids_clr <- do.call(rbind, lapply(paste0("C", 1:4), function(cl) {
  samples <- cluster$sample[cluster$cluster == cl]
  colMeans(tcga_clr[samples, , drop = FALSE], na.rm = TRUE)
}))
rownames(centroids_clr) <- paste0("C", 1:4)

scores <- data.frame(sample = rownames(external_fraction), stringsAsFactors = FALSE)
for (cl in paste0("C", 1:4)) {
  scores[[paste0(cl, "_spearman_raw")]] <- apply(external_fraction, 1, safe_cor, y = centroids[cl, ], method = "spearman")
  scores[[paste0(cl, "_pearson_asin")]] <- apply(external_asin, 1, safe_cor, y = centroids_asin[cl, ], method = "pearson")
  scores[[paste0(cl, "_pearson_clr")]] <- apply(external_clr, 1, safe_cor, y = centroids_clr[cl, ], method = "pearson")
  scores[[paste0(cl, "_cosine_raw")]] <- apply(external_fraction, 1, cosine_sim, y = centroids[cl, ])
}

scores$C4_minus_C3_spearman_raw <- scores$C4_spearman_raw - scores$C3_spearman_raw
scores$C4_minus_C3_pearson_asin <- scores$C4_pearson_asin - scores$C3_pearson_asin
scores$C4_minus_C3_pearson_clr <- scores$C4_pearson_clr - scores$C3_pearson_clr
scores$C4_minus_C3_cosine_raw <- scores$C4_cosine_raw - scores$C3_cosine_raw
scores$C4_minus_C2_spearman_raw <- scores$C4_spearman_raw - scores$C2_spearman_raw
scores$C2_minus_C3_spearman_raw <- scores$C2_spearman_raw - scores$C3_spearman_raw

for (n in c(5, 8, 10, 15, 20)) {
  c4_cells <- c4_high[seq_len(n)]
  c3_cells <- c3_high[seq_len(n)]
  scores[[paste0("auto_top", n, "_C4mean_minus_C3mean_raw")]] <-
    row_mean_existing(external_fraction, c4_cells) - row_mean_existing(external_fraction, c3_cells)
  scores[[paste0("auto_top", n, "_C4sum_over_C3sum_logratio")]] <-
    log((row_sum_existing(external_fraction, c4_cells) + 1e-5) /
          (row_sum_existing(external_fraction, c3_cells) + 1e-5))
  scores[[paste0("auto_top", n, "_C4mean_minus_C3mean_asin")]] <-
    row_mean_existing(external_asin, c4_cells) - row_mean_existing(external_asin, c3_cells)
}

effector_cells <- c(
  "CD8Tem", "CD8Temra", "CD8Trm", "Tc", "Tex", "CD8Tcm", "CD8Tn",
  "cNK", "NKreg", "NKT", "gdT", "MAIT",
  "Th1", "Tfh", "CD4Tem", "CD4Temra", "CD4Trm",
  "cDC1", "cDC2", "pDC", "PB", "PC", "MBC", "BGC"
)
suppressive_c3_cells <- c("MDSC", "M1", "M2", "TAM", "cMo", "intMo", "ncMo", "Neutrophil", "Treg")
scores$manual_effector_minus_C3suppressive_raw <-
  row_sum_existing(external_fraction, effector_cells) - row_sum_existing(external_fraction, suppressive_c3_cells)
scores$manual_effector_over_C3suppressive_logratio <-
  log((row_sum_existing(external_fraction, effector_cells) + 1e-5) /
        (row_sum_existing(external_fraction, suppressive_c3_cells) + 1e-5))
scores$manual_effector_mean_minus_C3suppressive_asin <-
  row_mean_existing(external_asin, effector_cells) - row_mean_existing(external_asin, suppressive_c3_cells)
scores$C3_MDSC_M1_negative <- -row_sum_existing(external_fraction, c("MDSC", "M1"))
scores$C4effector_over_MDSC_M1_logratio <-
  log((row_sum_existing(external_fraction, effector_cells) + 1e-5) /
        (row_sum_existing(external_fraction, c("MDSC", "M1")) + 1e-5))

# Train a C3-vs-C4 discriminant on TCGA only, then apply to the ICB samples.
train_idx <- cluster$sample[cluster$cluster %in% c("C3", "C4")]
train_y <- ifelse(cluster$cluster[match(train_idx, cluster$sample)] == "C4", 1L, 0L)
names(train_y) <- train_idx

for (tr_name in c("asin", "clr")) {
  x_train <- if (tr_name == "asin") tcga_asin[train_idx, , drop = FALSE] else tcga_clr[train_idx, , drop = FALSE]
  x_test <- if (tr_name == "asin") external_asin else external_clr
  xs <- z_by_train(x_train, x_test)
  fit_ridge <- suppressWarnings(cv.glmnet(xs$train, train_y, family = "binomial", alpha = 0, nfolds = 5))
  fit_lasso <- suppressWarnings(cv.glmnet(xs$train, train_y, family = "binomial", alpha = 1, nfolds = 5))
  scores[[paste0("TCGA_C4_vs_C3_glmnet_ridge_probability_", tr_name)]] <-
    as.numeric(predict(fit_ridge, newx = xs$test, s = "lambda.min", type = "response"))
  scores[[paste0("TCGA_C4_vs_C3_glmnet_lasso_probability_", tr_name)]] <-
    as.numeric(predict(fit_lasso, newx = xs$test, s = "lambda.min", type = "response"))

  fit_rf <- ranger(
    x = as.data.frame(xs$train, check.names = FALSE),
    y = factor(train_y, levels = c(0, 1)),
    probability = TRUE,
    classification = TRUE,
    num.trees = 1000,
    mtry = 15,
    min.node.size = 10,
    seed = 20260622,
    num.threads = min(8L, parallel::detectCores(logical = TRUE))
  )
  scores[[paste0("TCGA_C4_vs_C3_ranger_probability_", tr_name)]] <-
    predict(fit_rf, data = as.data.frame(xs$test, check.names = FALSE))$predictions[, "1"]
}

fwrite(scores, file.path(out_dir, "external_ICB_C3_C4_axis_scores.txt"), sep = "\t")

overall <- do.call(rbind, lapply(setdiff(colnames(scores), "sample"), function(sc) {
  cbind(score = sc, score_summary(y, scores[[sc]]))
}))
overall <- overall[order(overall$OrientedAUC, decreasing = TRUE), ]
fwrite(overall, file.path(out_dir, "overall_C3_C4_axis_score_auc.txt"), sep = "\t")

meta <- score_df
meta$Run <- as.character(meta$Run)
rownames(meta) <- meta$Run
subgroup_vars <- intersect(c("SRA_study", "Cancer", "Cancer_type", "Drug", "Anti_target", "disease", "Biopsy_Time"), colnames(meta))

subgroup_rows <- list()
for (sc in setdiff(colnames(scores), "sample")) {
  for (var in subgroup_vars) {
    vals <- as.character(meta[scores$sample, var])
    vals[is.na(vals) | vals == ""] <- "Unknown"
    for (lv in sort(unique(vals))) {
      idx <- which(vals == lv)
      yy <- y[idx]
      if (length(idx) < 20 || sum(yy == 1) < 5 || sum(yy == 0) < 5) next
      sm <- score_summary(yy, scores[[sc]][idx])
      subgroup_rows[[length(subgroup_rows) + 1L]] <- cbind(
        score = sc,
        subgroup_variable = var,
        subgroup = lv,
        N = length(idx),
        N_R = sum(yy == 1),
        N_NR = sum(yy == 0),
        sm,
        stringsAsFactors = FALSE
      )
    }
  }
}
subgroup_auc <- do.call(rbind, subgroup_rows)
subgroup_auc <- subgroup_auc[order(subgroup_auc$OrientedAUC, decreasing = TRUE), ]
fwrite(subgroup_auc, file.path(out_dir, "subgroup_C3_C4_axis_score_auc_min20_min5R5NR.txt"), sep = "\t")

study_rows <- subgroup_auc[subgroup_auc$subgroup_variable == "SRA_study", , drop = FALSE]
fwrite(study_rows, file.path(out_dir, "SRA_study_C3_C4_axis_score_auc_min20_min5R5NR.txt"), sep = "\t")

p1 <- ggplot(head(overall, 25), aes(x = reorder(score, OrientedAUC), y = OrientedAUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Oriented AUC", title = "Overall ICB response prediction by C3/C4-axis scores") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "overall_top25_C3_C4_axis_score_auc.pdf"), p1, width = 9, height = 6)
ggsave(file.path(out_dir, "overall_top25_C3_C4_axis_score_auc.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(head(study_rows, 30), aes(x = reorder(paste0(subgroup, " | ", score), OrientedAUC), y = OrientedAUC)) +
  geom_col(fill = "#F66463", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, color = "black", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = "Oriented AUC", title = "Best SRA-study subgroups for C3/C4-axis scores") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "SRA_study_top30_C3_C4_axis_score_auc.pdf"), p2, width = 10, height = 7)
ggsave(file.path(out_dir, "SRA_study_top30_C3_C4_axis_score_auc.png"), p2, width = 10, height = 7, dpi = 300)

writeLines(c(
  "C3/C4-axis scoring for ICB response",
  paste0("External samples: ", length(y), "; R=", sum(y == 1), "; NR=", sum(y == 0)),
  "Only ImmuCellAI2 cell-fraction features were used.",
  "No IFN-gamma, cytotoxicity, checkpoint, MHC, T cell inflamed GEP, cancer type, drug, target, study, or batch variables were used as model predictors.",
  "Subgroup AUC tables are exploratory external-validation screens; high AUC subgroups should be presented as hypothesis-generating unless pre-specified."
), file.path(out_dir, "README_C3_C4_axis_ICB_response_scores.txt"))

print(cluster_sizes)
print(head(top_table, 15))
print(head(overall, 20))
print(head(study_rows, 20))
print(head(subgroup_auc, 30))
