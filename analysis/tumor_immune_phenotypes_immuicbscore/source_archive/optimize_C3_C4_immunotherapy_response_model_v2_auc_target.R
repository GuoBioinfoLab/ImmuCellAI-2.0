options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
out_dir <- file.path(base_dir, "optimized_models_auc_target")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(nnet)
  library(gbm)
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

safe01 <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  x
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

remove_bad_columns <- function(x_train, x_test) {
  s <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  if (!any(keep)) {
    return(list(train = matrix(0, nrow(x_train), 1, dimnames = list(rownames(x_train), "zero")),
                test = matrix(0, nrow(x_test), 1, dimnames = list(rownames(x_test), "zero"))))
  }
  list(train = x_train[, keep, drop = FALSE], test = x_test[, keep, drop = FALSE])
}

scale_by_train <- function(x_train, x_test) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  center <- colMeans(x_train, na.rm = TRUE)
  scalev <- apply(x_train, 2, sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scalev, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scalev, "/")
  )
}

one_hot_meta <- function(meta, cols) {
  dat <- meta[, cols, drop = FALSE]
  for (cc in names(dat)) {
    dat[[cc]] <- as.character(dat[[cc]])
    dat[[cc]][is.na(dat[[cc]]) | dat[[cc]] == ""] <- "Unknown"
    dat[[cc]] <- factor(dat[[cc]])
  }
  mm <- model.matrix(~ . - 1, data = dat)
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  mm
}

build_feature_matrix <- function(feature_set, cell_mat, score_df) {
  cells <- safe01(cell_mat)
  cells_asin <- asin(sqrt(cells))
  colnames(cells_asin) <- paste0("cell_asin_", make.names(colnames(cells_asin)))

  cells_logit <- log((cells + 1e-5) / (1 - pmin(cells, 1 - 1e-5)))
  colnames(cells_logit) <- paste0("cell_logit_", make.names(colnames(cells_logit)))

  score_cols <- c(
    "C3_like_spearman", "C4_like_spearman",
    "C3_like_pearson", "C4_like_pearson",
    "C3_like_cosine", "C4_like_cosine",
    "Weighted_C4_minus_C3",
    "C4_minus_C3_spearman", "C4_minus_C3_pearson", "C4_minus_C3_cosine"
  )
  score_cols <- intersect(score_cols, colnames(score_df))
  scores <- as.matrix(score_df[, score_cols, drop = FALSE])
  storage.mode(scores) <- "numeric"
  scores[!is.finite(scores)] <- 0
  colnames(scores) <- paste0("score_", make.names(colnames(scores)))

  basic_meta_cols <- intersect(c("Cancer", "Cancer_type", "Tissue", "Drug", "Anti_target", "Biopsy_Time", "Gender", "disease"), colnames(score_df))
  study_meta_cols <- intersect(c(basic_meta_cols, "SRA_study", "SRA_Study", "batch"), colnames(score_df))
  basic_meta <- if (length(basic_meta_cols)) one_hot_meta(score_df, basic_meta_cols) else NULL
  study_meta <- if (length(study_meta_cols)) one_hot_meta(score_df, study_meta_cols) else NULL

  if (feature_set == "baseline_scores") return(scores)
  if (feature_set == "cells_asin") return(cells_asin)
  if (feature_set == "cells_logit") return(cells_logit)
  if (feature_set == "scores_cells_asin") return(cbind(scores, cells_asin))
  if (feature_set == "scores_cells_logit") return(cbind(scores, cells_logit))
  if (feature_set == "cells_asin_basic_meta") return(cbind(cells_asin, basic_meta))
  if (feature_set == "scores_cells_asin_basic_meta") return(cbind(scores, cells_asin, basic_meta))
  if (feature_set == "scores_cells_asin_study_meta") return(cbind(scores, cells_asin, study_meta))
  stop("Unknown feature_set: ", feature_set)
}

predict_glmnet <- function(y_train, x_train, x_test, alpha = 0.5) {
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(x_train, y_train, family = "binomial", alpha = alpha,
                               type.measure = "auc", nfolds = 5)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

predict_rf <- function(y_train, x_train, x_test, ntree = 500) {
  fit <- tryCatch(
    randomForest(x = as.data.frame(x_train), y = factor(y_train, levels = c(0, 1)),
                 ntree = ntree, importance = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newdata = as.data.frame(x_test), type = "prob")[, "1"])
}

predict_svm_radial <- function(y_train, x_train, x_test) {
  gamma_value <- 1 / max(1, ncol(x_train))
  fit <- tryCatch(
    svm(x = x_train, y = factor(y_train, levels = c(0, 1)), kernel = "radial",
        probability = TRUE, cost = 1, gamma = gamma_value, scale = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pred <- predict(fit, newdata = x_test, probability = TRUE)
  prob <- attr(pred, "probabilities")
  if (is.null(prob) || !"1" %in% colnames(prob)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(prob[, "1"])
}

predict_nnet <- function(y_train, x_train, x_test, size = 5, decay = 0.1) {
  dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train))
  fit <- tryCatch(
    suppressWarnings(nnet(y ~ ., data = dat, size = size, decay = decay,
                          maxit = 500, trace = FALSE, MaxNWts = 20000)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pr <- tryCatch(predict(fit, newdata = as.data.frame(x_test), type = "raw"),
                 error = function(e) NULL)
  if (is.null(pr)) return(rep(NA_real_, nrow(x_test)))
  if (is.matrix(pr) && "1" %in% colnames(pr)) return(as.numeric(pr[, "1"]))
  as.numeric(pr)
}

predict_gbm <- function(y_train, x_train, x_test, n_trees = 1200, depth = 2, shrinkage = 0.01) {
  dat <- data.frame(y = y_train, as.data.frame(x_train))
  fit <- tryCatch(
    suppressWarnings(gbm(y ~ ., data = dat, distribution = "bernoulli",
                         n.trees = n_trees, interaction.depth = depth,
                         shrinkage = shrinkage, bag.fraction = 0.75,
                         train.fraction = 1, verbose = FALSE)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newdata = as.data.frame(x_test), n.trees = n_trees, type = "response"))
}

predict_pca_glmnet <- function(y_train, x_train, x_test, n_pc = 10, alpha = 0.5) {
  n_pc <- min(n_pc, ncol(x_train), nrow(x_train) - 2L)
  if (n_pc < 1) return(rep(NA_real_, nrow(x_test)))
  pca <- tryCatch(prcomp(x_train, center = FALSE, scale. = FALSE), error = function(e) NULL)
  if (is.null(pca)) return(rep(NA_real_, nrow(x_test)))
  xt <- pca$x[, seq_len(n_pc), drop = FALSE]
  xv <- x_test %*% pca$rotation[, seq_len(n_pc), drop = FALSE]
  predict_glmnet(y_train, xt, xv, alpha = alpha)
}

fit_predict_model <- function(model_name, y_train, x_train, x_test) {
  if (model_name == "glmnet_ridge") return(predict_glmnet(y_train, x_train, x_test, alpha = 0))
  if (model_name == "glmnet_elastic") return(predict_glmnet(y_train, x_train, x_test, alpha = 0.5))
  if (model_name == "glmnet_lasso") return(predict_glmnet(y_train, x_train, x_test, alpha = 1))
  if (model_name == "random_forest") return(predict_rf(y_train, x_train, x_test, ntree = 500))
  if (model_name == "svm_radial") return(predict_svm_radial(y_train, x_train, x_test))
  if (model_name == "nnet5") return(predict_nnet(y_train, x_train, x_test, size = 5, decay = 0.1))
  if (model_name == "gbm_depth2") return(predict_gbm(y_train, x_train, x_test, n_trees = 1200, depth = 2, shrinkage = 0.01))
  if (model_name == "gbm_depth3") return(predict_gbm(y_train, x_train, x_test, n_trees = 1500, depth = 3, shrinkage = 0.008))
  if (model_name == "pca10_glmnet_elastic") return(predict_pca_glmnet(y_train, x_train, x_test, n_pc = 10, alpha = 0.5))
  if (model_name == "pca20_glmnet_elastic") return(predict_pca_glmnet(y_train, x_train, x_test, n_pc = 20, alpha = 0.5))
  stop("Unknown model_name: ", model_name)
}

run_repeated_cv <- function(x, y, model_name, label, n_repeats = 20, k = 5) {
  all_pred <- vector("list", n_repeats)
  auc_rows <- vector("list", n_repeats)
  for (rr in seq_len(n_repeats)) {
    folds <- make_stratified_folds(y, k = k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      train_idx <- which(folds != ff)
      test_idx <- which(folds == ff)
      xb <- remove_bad_columns(x[train_idx, , drop = FALSE], x[test_idx, , drop = FALSE])
      xs <- scale_by_train(xb$train, xb$test)
      pred[test_idx] <- fit_predict_model(model_name, y[train_idx], xs$train, xs$test)
    }
    all_pred[[rr]] <- data.frame(repeat_id = rr, sample = names(y), y = y, prediction = pred,
                                 model = label, stringsAsFactors = FALSE)
    auc_rows[[rr]] <- data.frame(repeat_id = rr, model = label,
                                 AUC = auc_manual(y, pred),
                                 Mean_R = mean(pred[y == 1], na.rm = TRUE),
                                 Mean_NR = mean(pred[y == 0], na.rm = TRUE),
                                 stringsAsFactors = FALSE)
    message(label, " repeat ", rr, "/", n_repeats, " AUC=", round(auc_rows[[rr]]$AUC, 4))
  }
  list(predictions = do.call(rbind, all_pred), auc = do.call(rbind, auc_rows))
}

message("Reading data...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
fraction_df <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
names(fraction_df)[1] <- "Run"
rownames(fraction_df) <- fraction_df$Run
fraction_df$Run <- NULL

response_binary <- rep(NA_integer_, nrow(score_df))
if ("ResponseGroup" %in% colnames(score_df)) {
  response_binary[score_df$ResponseGroup == "Responder"] <- 1L
  response_binary[score_df$ResponseGroup == "Non-responder"] <- 0L
}
if ("Response" %in% colnames(score_df)) {
  response_binary[is.na(response_binary) & score_df$Response == "R"] <- 1L
  response_binary[is.na(response_binary) & score_df$Response == "NR"] <- 0L
}
if ("ResponseBinary" %in% colnames(score_df)) {
  rb <- as.character(score_df$ResponseBinary)
  response_binary[is.na(response_binary) & rb == "1"] <- 1L
  response_binary[is.na(response_binary) & rb == "0"] <- 0L
}
score_df$ResponseBinaryClean <- response_binary
score_df <- score_df[!is.na(score_df$ResponseBinaryClean), , drop = FALSE]
score_df$Run <- as.character(score_df$Run)
score_df <- score_df[score_df$Run %in% rownames(fraction_df), , drop = FALSE]
cell_mat <- as.matrix(fraction_df[score_df$Run, , drop = FALSE])
y <- as.integer(score_df$ResponseBinaryClean)
names(y) <- score_df$Run

message("Samples with response: ", length(y), "; responders: ", sum(y == 1), "; nonresponders: ", sum(y == 0))

feature_sets <- c(
  "baseline_scores",
  "cells_asin",
  "scores_cells_asin",
  "scores_cells_logit",
  "cells_asin_basic_meta",
  "scores_cells_asin_basic_meta",
  "scores_cells_asin_study_meta"
)

model_grid <- list(
  strict = expand.grid(
    feature_set = c("baseline_scores", "scores_cells_asin", "scores_cells_logit"),
    model_name = c("glmnet_ridge", "glmnet_elastic", "pca20_glmnet_elastic", "random_forest", "gbm_depth2", "gbm_depth3"),
    stringsAsFactors = FALSE
  ),
  clinical_context = expand.grid(
    feature_set = c("cells_asin_basic_meta", "scores_cells_asin_basic_meta"),
    model_name = c("glmnet_ridge", "glmnet_elastic", "random_forest", "nnet5", "gbm_depth2", "gbm_depth3"),
    stringsAsFactors = FALSE
  ),
  study_adjusted_upper_bound = expand.grid(
    feature_set = "scores_cells_asin_study_meta",
    model_name = c("glmnet_ridge", "glmnet_elastic", "glmnet_lasso", "random_forest", "gbm_depth2"),
    stringsAsFactors = FALSE
  )
)

grid_df <- do.call(rbind, Map(function(mode, df) {
  df$mode <- mode
  df
}, names(model_grid), model_grid))

all_auc <- list()
all_predictions <- list()

for (ii in seq_len(nrow(grid_df))) {
  fs <- grid_df$feature_set[ii]
  mn <- grid_df$model_name[ii]
  mode <- grid_df$mode[ii]
  label <- paste(mode, fs, mn, sep = "__")
  message("Running ", label)
  x <- build_feature_matrix(fs, cell_mat, score_df)
  n_rep <- if (mode == "study_adjusted_upper_bound") 5 else 5
  res <- run_repeated_cv(x, y, mn, label, n_repeats = n_rep, k = 5)
  all_auc[[label]] <- res$auc
  all_predictions[[label]] <- res$predictions
}

auc_by_repeat <- do.call(rbind, all_auc)
predictions <- do.call(rbind, all_predictions)

summary_df <- aggregate(cbind(AUC, Mean_R, Mean_NR) ~ model, data = auc_by_repeat, FUN = mean)
summary_df$MedianAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = median)$AUC
summary_df$SDAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = sd)$AUC
summary_df$MaxAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = max)$AUC
summary_df$MinAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = min)$AUC
summary_df <- summary_df[order(summary_df$AUC, decreasing = TRUE), ]

fwrite(auc_by_repeat, file.path(out_dir, "auc_target_repeated5fold_auc_by_repeat.txt"), sep = "\t")
fwrite(predictions, file.path(out_dir, "auc_target_repeated5fold_predictions.txt"), sep = "\t")
fwrite(summary_df, file.path(out_dir, "auc_target_repeated5fold_auc_summary.txt"), sep = "\t")

best_model <- summary_df$model[1]
message("Best repeated-CV model: ", best_model, "; Mean AUC=", round(summary_df$AUC[1], 4))

best_pred <- predictions[predictions$model == best_model & predictions$repeat_id == 1, , drop = FALSE]
plot_df <- best_pred
plot_df$response <- factor(ifelse(plot_df$y == 1, "Responder", "Non-responder"),
                           levels = c("Non-responder", "Responder"))

p1 <- ggplot(auc_by_repeat, aes(x = reorder(model, AUC, FUN = median), y = AUC)) +
  geom_boxplot(fill = "#6BA5C3", color = "black", outlier.size = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#F66463") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Repeated 5-fold CV AUC",
       title = "Model optimization for immunotherapy response") +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "auc_target_repeated5fold_auc_boxplot.pdf"), p1, width = 10, height = 9)
ggsave(file.path(out_dir, "auc_target_repeated5fold_auc_boxplot.png"), p1, width = 10, height = 9, dpi = 300)

p2 <- ggplot(plot_df, aes(x = response, y = prediction, fill = response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
  scale_fill_manual(values = c("Non-responder" = "#95A8AC", "Responder" = "#F66463")) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Predicted response probability",
       title = paste0("Best model repeat 1: ", best_model),
       subtitle = paste0("AUC = ", round(auc_manual(plot_df$y, plot_df$prediction), 3))) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "auc_target_best_model_response_boxplot_repeat1.pdf"), p2, width = 4.2, height = 4.2)
ggsave(file.path(out_dir, "auc_target_best_model_response_boxplot_repeat1.png"), p2, width = 4.2, height = 4.2, dpi = 300)

writeLines(c(
  paste0("Best repeated 5-fold CV model: ", best_model),
  paste0("Mean AUC: ", signif(summary_df$AUC[1], 5)),
  paste0("Median AUC: ", signif(summary_df$MedianAUC[1], 5)),
  paste0("Max repeat AUC: ", signif(summary_df$MaxAUC[1], 5)),
  "",
  "Interpretation guide:",
  "strict = ImmuCellAI2 fractions and C3/C4 scores only.",
  "clinical_context = strict features plus cancer/drug/target/biopsy-time context.",
  "study_adjusted_upper_bound = additionally includes study/batch labels; treat as an optimistic upper bound because it can capture cohort effects."
), con = file.path(out_dir, "auc_target_best_model_readme.txt"))

print(summary_df)
