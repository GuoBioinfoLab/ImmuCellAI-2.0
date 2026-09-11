options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
out_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation", "optimized_models")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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

orient_score <- function(y_train, train_score, test_score) {
  a <- auc_manual(y_train, train_score)
  if (is.na(a)) return(test_score)
  if (a < 0.5) -test_score else test_score
}

make_stratified_folds <- function(y, k = 5L) {
  folds <- integer(length(y))
  for (yy in sort(unique(y))) {
    idx <- which(y == yy)
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

safe_scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train, na.rm = TRUE)
  scalev <- apply(x_train, 2, sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scalev, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scalev, "/"),
    center = center,
    scale = scalev
  )
}

fit_predict_glm <- function(y_train, x_train, x_test) {
  dat_train <- data.frame(y = y_train, x_train, check.names = FALSE)
  dat_test <- data.frame(x_test, check.names = FALSE)
  fit <- tryCatch(
    suppressWarnings(glm(y ~ ., data = dat_train, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pred <- tryCatch(
    suppressWarnings(predict(fit, newdata = dat_test, type = "response")),
    error = function(e) rep(NA_real_, nrow(x_test))
  )
  as.numeric(pred)
}

select_top_features <- function(y, x_train, top_n = 10L) {
  aucs <- apply(x_train, 2, function(v) auc_manual(y, v))
  score <- abs(aucs - 0.5)
  score[!is.finite(score)] <- -Inf
  names(sort(score, decreasing = TRUE))[seq_len(min(top_n, sum(is.finite(score))))]
}

predict_top_glm <- function(y_train, x_train, x_test, top_n = 10L, add_scores = NULL, add_scores_test = NULL) {
  features <- select_top_features(y_train, x_train, top_n = top_n)
  if (!length(features)) return(rep(NA_real_, nrow(x_test)))
  xt <- x_train[, features, drop = FALSE]
  xv <- x_test[, features, drop = FALSE]
  if (!is.null(add_scores)) {
    xt <- cbind(xt, add_scores)
    xv <- cbind(xv, add_scores_test)
  }
  fit_predict_glm(y_train, xt, xv)
}

predict_pca_glm <- function(y_train, x_train, x_test, n_pc = 10L, add_scores = NULL, add_scores_test = NULL) {
  n_pc <- min(n_pc, ncol(x_train), nrow(x_train) - 2L)
  if (n_pc < 1L) return(rep(NA_real_, nrow(x_test)))
  pca <- tryCatch(prcomp(x_train, center = FALSE, scale. = FALSE), error = function(e) NULL)
  if (is.null(pca)) return(rep(NA_real_, nrow(x_test)))
  xt <- pca$x[, seq_len(n_pc), drop = FALSE]
  xv <- x_test %*% pca$rotation[, seq_len(n_pc), drop = FALSE]
  colnames(xt) <- paste0("PC", seq_len(n_pc))
  colnames(xv) <- paste0("PC", seq_len(n_pc))
  if (!is.null(add_scores)) {
    xt <- cbind(xt, add_scores)
    xv <- cbind(xv, add_scores_test)
  }
  fit_predict_glm(y_train, as.data.frame(xt), as.data.frame(xv))
}

metric_from_prediction <- function(y, pred) {
  data.frame(
    AUC = auc_manual(y, pred),
    Mean_R = mean(pred[y == 1], na.rm = TRUE),
    Mean_NR = mean(pred[y == 0], na.rm = TRUE),
    Median_R = median(pred[y == 1], na.rm = TRUE),
    Median_NR = median(pred[y == 0], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

message("Reading scores and fraction matrix...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
frac <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
rownames(frac) <- frac[[1L]]
frac[[1L]] <- NULL
frac <- as.matrix(frac)
storage.mode(frac) <- "double"

dat <- score_df %>%
  filter(!is.na(ResponseBinary), Run %in% rownames(frac)) %>%
  mutate(ResponseBinary = as.integer(ResponseBinary))
frac <- frac[dat$Run, , drop = FALSE]

cell_features <- colnames(frac)
score_features <- c(
  "C4_minus_C3_spearman", "C3_like_spearman", "C4_like_spearman",
  "C4_minus_C3_pearson", "C4_minus_C3_cosine", "Weighted_C4_minus_C3"
)
score_mat <- as.matrix(dat[, score_features])
storage.mode(score_mat) <- "double"

# Compositional proportions: arcsine sqrt is a stable transform for fractions.
x_cell_raw <- asin(sqrt(pmax(pmin(frac, 1), 0)))
rownames(x_cell_raw) <- dat$Run

y <- dat$ResponseBinary
n <- length(y)
message("Samples with response: ", n, "; responders=", sum(y == 1), "; non-responders=", sum(y == 0))

models <- c(
  "baseline_C4_minus_C3_spearman",
  "baseline_negative_C3_like",
  "top5_cells_glm",
  "top10_cells_glm",
  "top20_cells_glm",
  "scores_plus_top10_cells_glm",
  "pca5_cells_glm",
  "pca10_cells_glm",
  "scores_plus_pca10_cells_glm",
  "all53_cells_glm",
  "scores_plus_all53_cells_glm"
)

n_repeats <- 30L
k_folds <- 5L
cv_predictions <- vector("list", n_repeats)
cv_summary <- vector("list", n_repeats)

for (rep_i in seq_len(n_repeats)) {
  message("CV repeat ", rep_i, "/", n_repeats)
  folds <- make_stratified_folds(y, k = k_folds)
  pred_mat <- matrix(NA_real_, nrow = n, ncol = length(models), dimnames = list(dat$Run, models))

  for (fold in seq_len(k_folds)) {
    test_idx <- which(folds == fold)
    train_idx <- setdiff(seq_len(n), test_idx)
    y_train <- y[train_idx]

    scaled <- safe_scale_train_test(x_cell_raw[train_idx, , drop = FALSE], x_cell_raw[test_idx, , drop = FALSE])
    x_train <- as.data.frame(scaled$train, check.names = FALSE)
    x_test <- as.data.frame(scaled$test, check.names = FALSE)

    score_scaled <- safe_scale_train_test(score_mat[train_idx, , drop = FALSE], score_mat[test_idx, , drop = FALSE])
    s_train <- as.data.frame(score_scaled$train, check.names = FALSE)
    s_test <- as.data.frame(score_scaled$test, check.names = FALSE)

    pred_mat[test_idx, "baseline_C4_minus_C3_spearman"] <-
      orient_score(y_train, score_mat[train_idx, "C4_minus_C3_spearman"], score_mat[test_idx, "C4_minus_C3_spearman"])
    pred_mat[test_idx, "baseline_negative_C3_like"] <-
      orient_score(y_train, -score_mat[train_idx, "C3_like_spearman"], -score_mat[test_idx, "C3_like_spearman"])

    pred_mat[test_idx, "top5_cells_glm"] <- predict_top_glm(y_train, x_train, x_test, top_n = 5)
    pred_mat[test_idx, "top10_cells_glm"] <- predict_top_glm(y_train, x_train, x_test, top_n = 10)
    pred_mat[test_idx, "top20_cells_glm"] <- predict_top_glm(y_train, x_train, x_test, top_n = 20)
    pred_mat[test_idx, "scores_plus_top10_cells_glm"] <- predict_top_glm(
      y_train, x_train, x_test, top_n = 10,
      add_scores = s_train[, c("C4_minus_C3_spearman", "C3_like_spearman"), drop = FALSE],
      add_scores_test = s_test[, c("C4_minus_C3_spearman", "C3_like_spearman"), drop = FALSE]
    )
    pred_mat[test_idx, "pca5_cells_glm"] <- predict_pca_glm(y_train, as.matrix(x_train), as.matrix(x_test), n_pc = 5)
    pred_mat[test_idx, "pca10_cells_glm"] <- predict_pca_glm(y_train, as.matrix(x_train), as.matrix(x_test), n_pc = 10)
    pred_mat[test_idx, "scores_plus_pca10_cells_glm"] <- predict_pca_glm(
      y_train, as.matrix(x_train), as.matrix(x_test), n_pc = 10,
      add_scores = s_train[, c("C4_minus_C3_spearman", "C3_like_spearman"), drop = FALSE],
      add_scores_test = s_test[, c("C4_minus_C3_spearman", "C3_like_spearman"), drop = FALSE]
    )
    pred_mat[test_idx, "all53_cells_glm"] <- fit_predict_glm(y_train, x_train, x_test)
    pred_mat[test_idx, "scores_plus_all53_cells_glm"] <- fit_predict_glm(y_train, cbind(x_train, s_train), cbind(x_test, s_test))
  }

  cv_predictions[[rep_i]] <- data.frame(Repeat = rep_i, Run = dat$Run, ResponseBinary = y, pred_mat, check.names = FALSE)
  cv_summary[[rep_i]] <- bind_rows(lapply(models, function(m) {
    cbind(Repeat = rep_i, Model = m, metric_from_prediction(y, pred_mat[, m]))
  }))
}

cv_pred_df <- bind_rows(cv_predictions)
cv_summary_df <- bind_rows(cv_summary)
cv_summary_stats <- cv_summary_df %>%
  group_by(Model) %>%
  summarise(
    MeanAUC = mean(AUC, na.rm = TRUE),
    MedianAUC = median(AUC, na.rm = TRUE),
    SDAUC = sd(AUC, na.rm = TRUE),
    Q25AUC = quantile(AUC, 0.25, na.rm = TRUE),
    Q75AUC = quantile(AUC, 0.75, na.rm = TRUE),
    Mean_R = mean(Mean_R, na.rm = TRUE),
    Mean_NR = mean(Mean_NR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(MeanAUC))

write.table(cv_pred_df, file.path(out_dir, "optimized_model_repeated5fold_predictions.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cv_summary_df, file.path(out_dir, "optimized_model_repeated5fold_auc_by_repeat.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cv_summary_stats, file.path(out_dir, "optimized_model_repeated5fold_auc_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Fitting final interpretable model on all response-labeled samples...")
scaled_all <- safe_scale_train_test(x_cell_raw, x_cell_raw)
x_all <- as.data.frame(scaled_all$train, check.names = FALSE)
top10 <- select_top_features(y, x_all, top_n = 10)
top20 <- select_top_features(y, x_all, top_n = 20)
write.table(data.frame(Top10 = top10), file.path(out_dir, "final_model_top10_selected_cells.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Top20 = top20), file.path(out_dir, "final_model_top20_selected_cells.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

final_train <- data.frame(y = y, x_all[, top10, drop = FALSE], check.names = FALSE)
final_fit <- suppressWarnings(glm(y ~ ., data = final_train, family = binomial()))
coef_df <- data.frame(
  Feature = names(coef(final_fit)),
  Coefficient = as.numeric(coef(final_fit)),
  stringsAsFactors = FALSE
)
write.table(coef_df, file.path(out_dir, "final_top10_cells_glm_coefficients.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Study-held-out validation for the best simple model family.
message("Running leave-study-out validation for top10 model...")
study_col <- if ("SRA_study" %in% colnames(dat)) "SRA_study" else "disease"
study_pred <- rep(NA_real_, n)
study_groups <- unique(dat[[study_col]])
for (g in study_groups) {
  test_idx <- which(dat[[study_col]] == g)
  train_idx <- setdiff(seq_len(n), test_idx)
  if (length(test_idx) < 5L || length(unique(y[test_idx])) < 2L || length(unique(y[train_idx])) < 2L) next
  scaled <- safe_scale_train_test(x_cell_raw[train_idx, , drop = FALSE], x_cell_raw[test_idx, , drop = FALSE])
  x_train <- as.data.frame(scaled$train, check.names = FALSE)
  x_test <- as.data.frame(scaled$test, check.names = FALSE)
  study_pred[test_idx] <- predict_top_glm(y[train_idx], x_train, x_test, top_n = 10)
}
study_auc <- auc_manual(y, study_pred)
study_summary <- data.frame(
  Model = "top10_cells_glm_leave_study_out",
  GroupColumn = study_col,
  N_predicted = sum(is.finite(study_pred)),
  AUC = study_auc,
  stringsAsFactors = FALSE
)
write.table(data.frame(Run = dat$Run, ResponseBinary = y, Study = dat[[study_col]], Prediction = study_pred),
            file.path(out_dir, "top10_cells_glm_leave_study_out_predictions.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(study_summary, file.path(out_dir, "top10_cells_glm_leave_study_out_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_auc <- ggplot(cv_summary_df, aes(x = reorder(Model, AUC, median), y = AUC, fill = Model)) +
  geom_boxplot(width = 0.62, outlier.size = 0.4, color = "black", linewidth = 0.25, show.legend = FALSE) +
  coord_flip() +
  theme_classic(base_size = 11) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)) +
  labs(x = NULL, y = "Repeated 5-fold CV AUC")

pdf(file.path(out_dir, "optimized_model_repeated5fold_auc_boxplot.pdf"), width = 8.2, height = 5.8, useDingbats = FALSE)
print(p_auc)
dev.off()
png(file.path(out_dir, "optimized_model_repeated5fold_auc_boxplot.png"), width = 2400, height = 1700, res = 300)
print(p_auc)
dev.off()

best_model <- cv_summary_stats$Model[1]
best_pred <- cv_pred_df %>%
  filter(Repeat == 1) %>%
  select(Run, ResponseBinary, all_of(best_model))
colnames(best_pred)[3] <- "Prediction"
p_best <- ggplot(best_pred, aes(x = factor(ResponseBinary, levels = c(0, 1), labels = c("Non-responder", "Responder")), y = Prediction,
                                fill = factor(ResponseBinary))) +
  geom_boxplot(width = 0.56, outlier.size = 0.6, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_jitter(width = 0.14, size = 0.45, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "#4E79A7FF", "1" = "#E15759FF")) +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)) +
  labs(x = NULL, y = paste0(best_model, " prediction"))

pdf(file.path(out_dir, "best_optimized_model_response_boxplot_repeat1.pdf"), width = 5.2, height = 4.6, useDingbats = FALSE)
print(p_best)
dev.off()
png(file.path(out_dir, "best_optimized_model_response_boxplot_repeat1.png"), width = 1560, height = 1380, res = 300)
print(p_best)
dev.off()

message("Done. Outputs written to: ", out_dir)
print(cv_summary_stats)
print(study_summary)
