options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
axis_dir <- file.path(base_dir, "C3_C4_axis_ICB_response_scores")
out_dir <- file.path(base_dir, "more_unified_ICB_score_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
axis_score_file <- file.path(axis_dir, "external_ICB_C3_C4_axis_scores.txt")
meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ranger)
  library(lightgbm)
  library(gbm)
  library(e1071)
  library(nnet)
  library(ggplot2)
})

set.seed(20260622)
n_threads <- min(8L, parallel::detectCores(logical = TRUE))

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

logit_transform <- function(x, eps = 1e-5) {
  p <- pmin(pmax(x, eps), 1 - eps)
  log(p / (1 - p))
}

row_sum_existing <- function(mat, cells) {
  cells <- intersect(cells, colnames(mat))
  if (!length(cells)) return(rep(NA_real_, nrow(mat)))
  rowSums(mat[, cells, drop = FALSE], na.rm = TRUE)
}

make_folds <- function(y, k = 5L) {
  folds <- integer(length(y))
  for (yy in sort(unique(y))) {
    idx <- sample(which(y == yy))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

scale_train_test <- function(x_train, x_test) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  x_train[!is.finite(x_train)] <- NA_real_
  x_test[!is.finite(x_test)] <- NA_real_
  med <- apply(x_train, 2, median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(x_train))) {
    x_train[!is.finite(x_train[, j]), j] <- med[j]
    x_test[!is.finite(x_test[, j]), j] <- med[j]
  }
  s <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
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

predict_model <- function(model_name, x_train, y_train, x_test) {
  if (model_name == "glmnet_ridge") {
    fit <- suppressWarnings(cv.glmnet(x_train, y_train, family = "binomial", alpha = 0,
                                      weights = class_weights(y_train), type.measure = "auc",
                                      nfolds = 5, standardize = FALSE))
    return(as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response")))
  }
  if (model_name == "glmnet_elastic") {
    fit <- suppressWarnings(cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5,
                                      weights = class_weights(y_train), type.measure = "auc",
                                      nfolds = 5, standardize = FALSE))
    return(as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response")))
  }
  if (model_name == "ranger_balanced") {
    fit <- ranger(
      x = as.data.frame(x_train, check.names = FALSE),
      y = factor(y_train, levels = c(0, 1)),
      probability = TRUE,
      classification = TRUE,
      num.trees = 1200,
      mtry = max(1, floor(sqrt(ncol(x_train)))),
      min.node.size = 8,
      class.weights = c("0" = class_weights(y_train)[which(y_train == 0)[1]],
                        "1" = class_weights(y_train)[which(y_train == 1)[1]]),
      seed = sample.int(1e6, 1),
      num.threads = n_threads
    )
    return(as.numeric(predict(fit, data = as.data.frame(x_test, check.names = FALSE))$predictions[, "1"]))
  }
  if (model_name == "ranger_wide") {
    fit <- ranger(
      x = as.data.frame(x_train, check.names = FALSE),
      y = factor(y_train, levels = c(0, 1)),
      probability = TRUE,
      classification = TRUE,
      num.trees = 1200,
      mtry = min(25, ncol(x_train)),
      min.node.size = 12,
      class.weights = c("0" = class_weights(y_train)[which(y_train == 0)[1]],
                        "1" = class_weights(y_train)[which(y_train == 1)[1]]),
      seed = sample.int(1e6, 1),
      num.threads = n_threads
    )
    return(as.numeric(predict(fit, data = as.data.frame(x_test, check.names = FALSE))$predictions[, "1"]))
  }
  if (model_name == "lightgbm_shallow") {
    dtrain <- lgb.Dataset(data = x_train, label = y_train, weight = class_weights(y_train))
    fit <- lgb.train(
      params = list(
        objective = "binary",
        metric = "auc",
        learning_rate = 0.03,
        num_leaves = 7,
        max_depth = 2,
        min_data_in_leaf = 25,
        feature_fraction = 0.85,
        bagging_fraction = 0.85,
        bagging_freq = 1,
        lambda_l1 = 0.05,
        lambda_l2 = 1,
        verbosity = -1,
        force_col_wise = TRUE,
        num_threads = n_threads
      ),
      data = dtrain,
      nrounds = 500,
      verbose = -1
    )
    return(as.numeric(predict(fit, x_test)))
  }
  if (model_name == "lightgbm_deeper") {
    dtrain <- lgb.Dataset(data = x_train, label = y_train, weight = class_weights(y_train))
    fit <- lgb.train(
      params = list(
        objective = "binary",
        metric = "auc",
        learning_rate = 0.02,
        num_leaves = 15,
        max_depth = 3,
        min_data_in_leaf = 20,
        feature_fraction = 0.85,
        bagging_fraction = 0.80,
        bagging_freq = 1,
        lambda_l1 = 0.1,
        lambda_l2 = 1,
        verbosity = -1,
        force_col_wise = TRUE,
        num_threads = n_threads
      ),
      data = dtrain,
      nrounds = 700,
      verbose = -1
    )
    return(as.numeric(predict(fit, x_test)))
  }
  if (model_name == "gbm_depth2") {
    dat <- data.frame(y = y_train, as.data.frame(x_train, check.names = FALSE))
    fit <- suppressWarnings(gbm(
      y ~ ., data = dat, distribution = "bernoulli",
      weights = class_weights(y_train),
      n.trees = 1200, interaction.depth = 2,
      shrinkage = 0.01, bag.fraction = 0.8,
      train.fraction = 1, verbose = FALSE
    ))
    return(as.numeric(predict(fit, newdata = as.data.frame(x_test, check.names = FALSE),
                              n.trees = 1200, type = "response")))
  }
  if (model_name == "svm_radial") {
    fit <- svm(x = x_train, y = factor(y_train, levels = c(0, 1)),
               kernel = "radial", probability = TRUE, cost = 1,
               gamma = 1 / max(1, ncol(x_train)), scale = FALSE,
               class.weights = c("0" = class_weights(y_train)[which(y_train == 0)[1]],
                                 "1" = class_weights(y_train)[which(y_train == 1)[1]]))
    pr <- predict(fit, newdata = x_test, probability = TRUE)
    prob <- attr(pr, "probabilities")
    return(as.numeric(prob[, "1"]))
  }
  if (model_name == "nnet4") {
    dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train, check.names = FALSE))
    fit <- suppressWarnings(nnet(y ~ ., data = dat, size = 4, decay = 0.1,
                                 weights = class_weights(y_train), maxit = 500,
                                 trace = FALSE, MaxNWts = 40000))
    pr <- predict(fit, newdata = as.data.frame(x_test, check.names = FALSE), type = "raw")
    if (is.matrix(pr) && "1" %in% colnames(pr)) return(as.numeric(pr[, "1"]))
    return(as.numeric(pr))
  }
  stop("Unknown model: ", model_name)
}

run_pooled_cv <- function(x, y, model_name, n_repeats = 5L, k = 5L) {
  rows <- vector("list", n_repeats)
  for (rr in seq_len(n_repeats)) {
    folds <- make_folds(y, k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      train <- which(folds != ff)
      test <- which(folds == ff)
      xs <- scale_train_test(x[train, , drop = FALSE], x[test, , drop = FALSE])
      pred[test] <- tryCatch(predict_model(model_name, xs$train, y[train], xs$test),
                             error = function(e) rep(NA_real_, length(test)))
    }
    rows[[rr]] <- data.frame(repeat_id = rr, AUC = auc_manual(y, pred),
                             Mean_R = mean(pred[y == 1], na.rm = TRUE),
                             Mean_NR = mean(pred[y == 0], na.rm = TRUE))
  }
  rbindlist(rows)
}

run_loso <- function(x, y, study, model_name, valid_studies) {
  rows <- list()
  preds <- list()
  for (st in valid_studies) {
    test <- which(study == st)
    train <- which(study != st)
    xs <- scale_train_test(x[train, , drop = FALSE], x[test, , drop = FALSE])
    pred <- tryCatch(predict_model(model_name, xs$train, y[train], xs$test),
                     error = function(e) rep(NA_real_, length(test)))
    rows[[length(rows) + 1L]] <- data.frame(
      SRA_study = st,
      N = length(test),
      N_R = sum(y[test] == 1),
      N_NR = sum(y[test] == 0),
      AUC = auc_manual(y[test], pred),
      Mean_R = mean(pred[y[test] == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y[test] == 0], na.rm = TRUE)
    )
    preds[[length(preds) + 1L]] <- data.frame(
      sample = names(y)[test],
      SRA_study = st,
      y = y[test],
      prediction = pred
    )
  }
  list(auc = rbindlist(rows), pred = rbindlist(preds))
}

message("Reading data...")
fraction <- read_fraction(fraction_file)
axis_scores <- fread(axis_score_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)

response <- rep(NA_integer_, nrow(meta))
response[meta$ResponseGroup == "Responder"] <- 1L
response[meta$ResponseGroup == "Non-responder"] <- 0L
response[is.na(response) & meta$Response == "R"] <- 1L
response[is.na(response) & meta$Response == "NR"] <- 0L
meta$response_binary <- response
meta <- meta[!is.na(meta$response_binary) & meta$Run %in% rownames(fraction) & meta$Run %in% axis_scores$sample, , drop = FALSE]
fraction <- fraction[meta$Run, , drop = FALSE]
axis_scores <- axis_scores[match(meta$Run, axis_scores$sample), , drop = FALSE]
y <- meta$response_binary
names(y) <- meta$Run
study <- as.character(meta$SRA_study)
message("Samples: ", length(y), "; R=", sum(y == 1), "; NR=", sum(y == 0))

raw <- fraction
asin <- asin(sqrt(pmax(raw, 0)))
clr <- clr_transform(raw)
logit <- logit_transform(raw)
colnames(asin) <- paste0("asin_", make.names(colnames(fraction), unique = TRUE))
colnames(clr) <- paste0("clr_", make.names(colnames(fraction), unique = TRUE))
colnames(logit) <- paste0("logit_", make.names(colnames(fraction), unique = TRUE))

axis_keep <- intersect(c(
  "C1_spearman_raw", "C2_spearman_raw", "C3_spearman_raw", "C4_spearman_raw",
  "C1_pearson_clr", "C2_pearson_clr", "C3_pearson_clr", "C4_pearson_clr",
  "C4_minus_C3_spearman_raw", "C4_minus_C3_pearson_clr", "C4_minus_C3_cosine_raw",
  "auto_top5_C4sum_over_C3sum_logratio", "auto_top10_C4sum_over_C3sum_logratio",
  "manual_effector_minus_C3suppressive_raw",
  "manual_effector_over_C3suppressive_logratio",
  "C3_MDSC_M1_negative",
  "C4effector_over_MDSC_M1_logratio",
  "TCGA_C4_vs_C3_ranger_probability_asin",
  "TCGA_C4_vs_C3_ranger_probability_clr"
), colnames(axis_scores))
axis_mat <- as.matrix(axis_scores[, axis_keep, drop = FALSE])
storage.mode(axis_mat) <- "numeric"
colnames(axis_mat) <- paste0("axis_", make.names(colnames(axis_mat), unique = TRUE))

effector_cells <- c("CD8Tem", "CD8Temra", "CD8Trm", "Tc", "Tex", "CD8Tcm", "CD8Tn",
                    "cNK", "NKreg", "NKT", "gdT", "MAIT", "Th1", "Tfh",
                    "CD4Tem", "CD4Temra", "CD4Trm", "cDC1", "cDC2", "pDC", "PB", "PC", "MBC", "BGC")
suppressive_cells <- c("MDSC", "M1", "M2", "TAM", "cMo", "intMo", "ncMo", "Neutrophil", "Treg")
myeloid_cells <- c("MDSC", "M1", "M2", "TAM", "cMo", "intMo", "ncMo", "Neutrophil", "M0")
tcell_cells <- c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg")
bio <- data.frame(
  effector_sum = row_sum_existing(fraction, effector_cells),
  suppressive_sum = row_sum_existing(fraction, suppressive_cells),
  myeloid_sum = row_sum_existing(fraction, myeloid_cells),
  tcell_sum = row_sum_existing(fraction, tcell_cells),
  MDSC = row_sum_existing(fraction, "MDSC"),
  M1 = row_sum_existing(fraction, "M1"),
  TAM = row_sum_existing(fraction, "TAM"),
  Treg = row_sum_existing(fraction, "Treg"),
  effector_over_suppressive = log((row_sum_existing(fraction, effector_cells) + 1e-5) /
                                    (row_sum_existing(fraction, suppressive_cells) + 1e-5)),
  tcell_over_myeloid = log((row_sum_existing(fraction, tcell_cells) + 1e-5) /
                             (row_sum_existing(fraction, myeloid_cells) + 1e-5))
)
colnames(bio) <- paste0("bio_", make.names(colnames(bio), unique = TRUE))
bio <- as.matrix(bio)

feature_sets <- list(
  axis_only = axis_mat,
  biology_compact = bio,
  cell_clr53 = clr,
  cell_asin53 = asin,
  axis_plus_clr53 = cbind(axis_mat, clr),
  axis_bio_clr53 = cbind(axis_mat, bio, clr),
  all_transforms_axis_bio = cbind(axis_mat, bio, asin, clr, logit)
)

valid_studies <- names(which(table(study) >= 20))
valid_studies <- valid_studies[sapply(valid_studies, function(st) {
  yy <- y[study == st]
  sum(yy == 1) >= 5 && sum(yy == 0) >= 5
})]

model_grid <- expand.grid(
  feature_set = names(feature_sets),
  model = c("glmnet_ridge", "glmnet_elastic", "ranger_balanced", "ranger_wide",
            "lightgbm_shallow", "lightgbm_deeper", "gbm_depth2", "svm_radial", "nnet4"),
  stringsAsFactors = FALSE
)

# Keep very slow/high-dimensional combinations out of the first broad run.
model_grid <- model_grid[!(model_grid$model %in% c("svm_radial", "nnet4") &
                             model_grid$feature_set %in% c("axis_bio_clr53", "all_transforms_axis_bio", "axis_plus_clr53")), ]

pooled_rows <- list()
loso_rows <- list()
pred_rows <- list()

pooled_file <- file.path(out_dir, "pooled_repeated5fold_auc_running.txt")
loso_file <- file.path(out_dir, "LOSO_auc_by_study_running.txt")
if (file.exists(pooled_file)) file.remove(pooled_file)
if (file.exists(loso_file)) file.remove(loso_file)

for (ii in seq_len(nrow(model_grid))) {
  fs <- model_grid$feature_set[ii]
  mn <- model_grid$model[ii]
  label <- paste(fs, mn, sep = "__")
  message("Running ", label)
  x <- as.matrix(feature_sets[[fs]])
  pooled <- run_pooled_cv(x, y, mn, n_repeats = 3, k = 5)
  pooled$model <- label
  pooled$feature_set <- fs
  pooled$method <- mn
  fwrite(pooled, pooled_file, sep = "\t", append = file.exists(pooled_file))
  pooled_rows[[label]] <- pooled

  loso <- run_loso(x, y, study, mn, valid_studies)
  loso$auc$model <- label
  loso$auc$feature_set <- fs
  loso$auc$method <- mn
  loso$pred$model <- label
  loso$pred$feature_set <- fs
  loso$pred$method <- mn
  fwrite(loso$auc, loso_file, sep = "\t", append = file.exists(loso_file))
  loso_rows[[label]] <- loso$auc
  pred_rows[[label]] <- loso$pred
}

pooled_all <- rbindlist(pooled_rows, fill = TRUE)
loso_all <- rbindlist(loso_rows, fill = TRUE)
pred_all <- rbindlist(pred_rows, fill = TRUE)

pooled_summary <- pooled_all[, .(
  PooledMeanAUC = mean(AUC, na.rm = TRUE),
  PooledMedianAUC = median(AUC, na.rm = TRUE),
  PooledSDAUC = sd(AUC, na.rm = TRUE),
  PooledMaxAUC = max(AUC, na.rm = TRUE)
), by = .(model, feature_set, method)]

loso_summary <- loso_all[, .(
  LOSO_MeanAUC = mean(AUC, na.rm = TRUE),
  LOSO_MedianAUC = median(AUC, na.rm = TRUE),
  LOSO_WeightedAUC = weighted.mean(AUC, N, na.rm = TRUE),
  LOSO_MinAUC = min(AUC, na.rm = TRUE),
  LOSO_MaxAUC = max(AUC, na.rm = TRUE),
  N_studies = .N,
  N_AUC_ge_0.65 = sum(AUC >= 0.65, na.rm = TRUE),
  N_AUC_ge_0.70 = sum(AUC >= 0.70, na.rm = TRUE),
  N_AUC_ge_0.75 = sum(AUC >= 0.75, na.rm = TRUE),
  N_AUC_ge_0.80 = sum(AUC >= 0.80, na.rm = TRUE)
), by = .(model, feature_set, method)]

summary_all <- merge(pooled_summary, loso_summary, by = c("model", "feature_set", "method"), all = TRUE)
summary_all <- summary_all[order(-summary_all$N_AUC_ge_0.75, -summary_all$LOSO_MeanAUC, -summary_all$PooledMeanAUC), ]

fwrite(pooled_all, file.path(out_dir, "pooled_repeated5fold_auc_by_repeat.txt"), sep = "\t")
fwrite(loso_all, file.path(out_dir, "LOSO_auc_by_study.txt"), sep = "\t")
fwrite(pred_all, file.path(out_dir, "LOSO_predictions.txt"), sep = "\t")
fwrite(summary_all, file.path(out_dir, "more_models_summary_ranked.txt"), sep = "\t")

best_loso <- summary_all$model[1]
best_overall <- summary_all$model[which.max(summary_all$PooledMeanAUC)]
fwrite(loso_all[loso_all$model == best_loso, ], file.path(out_dir, "BEST_by_LOSO_auc_by_study.txt"), sep = "\t")
fwrite(loso_all[loso_all$model == best_overall, ], file.path(out_dir, "BEST_by_pooled_auc_LOSO_by_study.txt"), sep = "\t")

p1 <- ggplot(head(summary_all, 25), aes(x = reorder(model, LOSO_MeanAUC), y = LOSO_MeanAUC, fill = N_AUC_ge_0.75)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  scale_fill_gradient(low = "#95A8AC", high = "#F66463") +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = "Mean LOSO AUC", fill = "N AUC >= 0.75",
       title = "Unified ICB score models ranked by LOSO performance")
ggsave(file.path(out_dir, "more_models_LOSO_summary_top25.pdf"), p1, width = 9, height = 6)
ggsave(file.path(out_dir, "more_models_LOSO_summary_top25.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(head(summary_all[order(-summary_all$PooledMeanAUC), ], 25),
             aes(x = reorder(model, PooledMeanAUC), y = PooledMeanAUC, fill = LOSO_MeanAUC)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  scale_fill_gradient(low = "#95A8AC", high = "#379DA5") +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = "Pooled repeated 5-fold CV AUC", fill = "Mean LOSO AUC",
       title = "Unified ICB score models ranked by pooled overall AUC")
ggsave(file.path(out_dir, "more_models_pooled_summary_top25.pdf"), p2, width = 9, height = 6)
ggsave(file.path(out_dir, "more_models_pooled_summary_top25.png"), p2, width = 9, height = 6, dpi = 300)

writeLines(c(
  paste0("Best model by LOSO AUC>=0.75 count: ", best_loso),
  paste0("Best model by pooled repeated-CV overall AUC: ", best_overall),
  "Predictors use only ImmuCellAI2 cell fractions and C3/C4 axis features derived from cell fractions.",
  "No IFN-gamma, cytotoxicity, checkpoint, MHC, T cell-inflamed GEP, cancer type, drug, target, study, or batch predictors were used.",
  "Pooled CV can overestimate cross-study performance; LOSO is stricter and should be prioritized for external validation claims."
), file.path(out_dir, "README_more_unified_ICB_score_models.txt"))

print(head(summary_all, 20))
cat("\nBest LOSO model by study:\n")
print(loso_all[loso_all$model == best_loso, ][order(-AUC), ])
cat("\nBest pooled-overall model by study:\n")
print(loso_all[loso_all$model == best_overall, ][order(-AUC), ])
