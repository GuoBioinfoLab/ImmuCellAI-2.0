options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
out_dir <- file.path(base_dir, "cell_fraction_only_model_optimization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(ranger)
  library(randomForest)
  library(gbm)
  library(nnet)
  library(e1071)
  library(kernlab)
  library(lightgbm)
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

make_stratified_folds <- function(y, k = 5L) {
  folds <- integer(length(y))
  for (yy in sort(unique(y))) {
    idx <- sample(which(y == yy))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
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

drop_constant <- function(x_train, x_test) {
  s <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  if (!any(keep)) {
    col <- matrix(0, nrow(x_train), 1, dimnames = list(rownames(x_train), "constant0"))
    col2 <- matrix(0, nrow(x_test), 1, dimnames = list(rownames(x_test), "constant0"))
    return(list(train = col, test = col2))
  }
  list(train = x_train[, keep, drop = FALSE], test = x_test[, keep, drop = FALSE])
}

select_top_features <- function(y, x_train, n = 20L) {
  aucs <- apply(x_train, 2, function(v) auc_manual(y, v))
  score <- abs(aucs - 0.5)
  score[!is.finite(score)] <- -Inf
  names(sort(score, decreasing = TRUE))[seq_len(min(n, sum(is.finite(score))))]
}

build_transform <- function(cell_mat, transform_name) {
  x <- as.matrix(cell_mat)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  eps <- 1e-5
  if (transform_name == "raw") {
    z <- x
  } else if (transform_name == "asin") {
    z <- asin(sqrt(x))
  } else if (transform_name == "sqrt") {
    z <- sqrt(x)
  } else if (transform_name == "logit") {
    p <- pmin(pmax(x, eps), 1 - eps)
    z <- log(p / (1 - p))
  } else if (transform_name == "clr") {
    lx <- log(x + eps)
    z <- sweep(lx, 1, rowMeans(lx), "-")
  } else if (transform_name == "rank") {
    z <- t(apply(x, 1, function(v) rank(v, ties.method = "average") / length(v)))
    colnames(z) <- colnames(x)
    rownames(z) <- rownames(x)
  } else {
    stop("Unknown transform: ", transform_name)
  }
  colnames(z) <- paste0(transform_name, "_", make.names(colnames(z)))
  z
}

predict_glmnet <- function(y_train, x_train, x_test, alpha = 0.5) {
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(x_train, y_train, family = "binomial",
                               alpha = alpha, type.measure = "auc",
                               nfolds = 5, parallel = FALSE)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

predict_ranger <- function(y_train, x_train, x_test, mtry = NULL, min_node = 5L) {
  dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train, check.names = FALSE))
  if (is.null(mtry)) mtry <- max(1L, floor(sqrt(ncol(x_train))))
  fit <- tryCatch(
    ranger(y ~ ., data = dat, probability = TRUE, num.trees = 1200,
           mtry = mtry, min.node.size = min_node,
           classification = TRUE, seed = sample.int(1e6, 1),
           num.threads = min(8L, parallel::detectCores(logical = TRUE))),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pr <- predict(fit, data = as.data.frame(x_test, check.names = FALSE))$predictions
  as.numeric(pr[, "1"])
}

predict_lgb <- function(y_train, x_train, x_test, leaves = 15L, depth = 3L, lr = 0.03, rounds = 500L) {
  dtrain <- tryCatch(lgb.Dataset(data = x_train, label = y_train), error = function(e) NULL)
  if (is.null(dtrain)) return(rep(NA_real_, nrow(x_test)))
  params <- list(
    objective = "binary",
    metric = "auc",
    learning_rate = lr,
    num_leaves = leaves,
    max_depth = depth,
    min_data_in_leaf = 20,
    feature_fraction = 0.85,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    lambda_l1 = 0.05,
    lambda_l2 = 1,
    verbosity = -1,
    num_threads = min(8L, parallel::detectCores(logical = TRUE)),
    force_col_wise = TRUE
  )
  fit <- tryCatch(
    lgb.train(params = params, data = dtrain, nrounds = rounds, verbose = -1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, x_test))
}

predict_svm <- function(y_train, x_train, x_test, cost = 1, gamma = NULL) {
  if (is.null(gamma)) gamma <- 1 / max(1, ncol(x_train))
  fit <- tryCatch(
    svm(x = x_train, y = factor(y_train, levels = c(0, 1)),
        kernel = "radial", probability = TRUE, cost = cost,
        gamma = gamma, scale = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pred <- predict(fit, newdata = x_test, probability = TRUE)
  prob <- attr(pred, "probabilities")
  if (is.null(prob) || !"1" %in% colnames(prob)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(prob[, "1"])
}

predict_gbm <- function(y_train, x_train, x_test, depth = 2L, trees = 1200L) {
  dat <- data.frame(y = y_train, as.data.frame(x_train, check.names = FALSE))
  fit <- tryCatch(
    suppressWarnings(gbm(y ~ ., data = dat, distribution = "bernoulli",
                         n.trees = trees, interaction.depth = depth,
                         shrinkage = 0.01, bag.fraction = 0.75,
                         train.fraction = 1, verbose = FALSE)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newdata = as.data.frame(x_test, check.names = FALSE),
                     n.trees = trees, type = "response"))
}

predict_nnet <- function(y_train, x_train, x_test, size = 4L, decay = 0.1) {
  dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train, check.names = FALSE))
  fit <- tryCatch(
    suppressWarnings(nnet(y ~ ., data = dat, size = size, decay = decay,
                          maxit = 500, trace = FALSE, MaxNWts = 30000)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pr <- tryCatch(predict(fit, newdata = as.data.frame(x_test, check.names = FALSE), type = "raw"),
                 error = function(e) NULL)
  if (is.null(pr)) return(rep(NA_real_, nrow(x_test)))
  if (is.matrix(pr) && "1" %in% colnames(pr)) return(as.numeric(pr[, "1"]))
  as.numeric(pr)
}

predict_model <- function(model_name, y_train, x_train, x_test) {
  if (model_name == "glmnet_ridge") return(predict_glmnet(y_train, x_train, x_test, alpha = 0))
  if (model_name == "glmnet_elastic") return(predict_glmnet(y_train, x_train, x_test, alpha = 0.5))
  if (model_name == "glmnet_lasso") return(predict_glmnet(y_train, x_train, x_test, alpha = 1))
  if (model_name == "ranger_sqrt_node5") return(predict_ranger(y_train, x_train, x_test, mtry = floor(sqrt(ncol(x_train))), min_node = 5))
  if (model_name == "ranger_mtry15_node5") return(predict_ranger(y_train, x_train, x_test, mtry = min(15, ncol(x_train)), min_node = 5))
  if (model_name == "ranger_mtry25_node10") return(predict_ranger(y_train, x_train, x_test, mtry = min(25, ncol(x_train)), min_node = 10))
  if (model_name == "lightgbm_d2_l7") return(predict_lgb(y_train, x_train, x_test, leaves = 7, depth = 2, lr = 0.03, rounds = 700))
  if (model_name == "lightgbm_d3_l15") return(predict_lgb(y_train, x_train, x_test, leaves = 15, depth = 3, lr = 0.03, rounds = 700))
  if (model_name == "lightgbm_d4_l31") return(predict_lgb(y_train, x_train, x_test, leaves = 31, depth = 4, lr = 0.02, rounds = 900))
  if (model_name == "svm_radial_c1") return(predict_svm(y_train, x_train, x_test, cost = 1))
  if (model_name == "svm_radial_c4") return(predict_svm(y_train, x_train, x_test, cost = 4))
  if (model_name == "gbm_d2") return(predict_gbm(y_train, x_train, x_test, depth = 2, trees = 1200))
  if (model_name == "gbm_d3") return(predict_gbm(y_train, x_train, x_test, depth = 3, trees = 1500))
  if (model_name == "nnet4") return(predict_nnet(y_train, x_train, x_test, size = 4, decay = 0.1))
  stop("Unknown model: ", model_name)
}

run_cv_model <- function(x_all, y, transform_name, model_name, top_n = NA_integer_,
                         n_repeats = 5L, k = 5L) {
  rows_auc <- vector("list", n_repeats)
  rows_pred <- vector("list", n_repeats)
  label <- paste(transform_name, model_name, ifelse(is.na(top_n), "all53", paste0("top", top_n)), sep = "__")
  for (rr in seq_len(n_repeats)) {
    folds <- make_stratified_folds(y, k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      train_idx <- which(folds != ff)
      test_idx <- which(folds == ff)
      xb <- drop_constant(x_all[train_idx, , drop = FALSE], x_all[test_idx, , drop = FALSE])
      x_train <- xb$train
      x_test <- xb$test
      if (!is.na(top_n)) {
        top_features <- select_top_features(y[train_idx], x_train, n = top_n)
        x_train <- x_train[, top_features, drop = FALSE]
        x_test <- x_test[, top_features, drop = FALSE]
      }
      xs <- scale_by_train(x_train, x_test)
      pred[test_idx] <- predict_model(model_name, y[train_idx], xs$train, xs$test)
    }
    auc <- auc_manual(y, pred)
    rows_auc[[rr]] <- data.frame(model = label, transform = transform_name,
                                 method = model_name, top_n = ifelse(is.na(top_n), "all53", paste0("top", top_n)),
                                 repeat_id = rr, AUC = auc,
                                 Mean_R = mean(pred[y == 1], na.rm = TRUE),
                                 Mean_NR = mean(pred[y == 0], na.rm = TRUE),
                                 stringsAsFactors = FALSE)
    rows_pred[[rr]] <- data.frame(model = label, repeat_id = rr, sample = names(y),
                                  y = y, prediction = pred, stringsAsFactors = FALSE)
    message(label, " repeat ", rr, "/", n_repeats, " AUC=", round(auc, 4))
  }
  list(auc = do.call(rbind, rows_auc), predictions = do.call(rbind, rows_pred))
}

message("Reading response and ImmuCellAI2 cell fractions...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
fraction_df <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
names(fraction_df)[1] <- "Run"
rownames(fraction_df) <- as.character(fraction_df$Run)
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
storage.mode(cell_mat) <- "numeric"
y <- as.integer(score_df$ResponseBinaryClean)
names(y) <- score_df$Run

message("Samples: ", length(y), "; responders: ", sum(y == 1), "; nonresponders: ", sum(y == 0),
        "; cell features: ", ncol(cell_mat))

transforms <- c("asin", "logit", "clr", "sqrt", "raw", "rank")
model_grid <- expand.grid(
  transform = transforms,
  method = c(
    "glmnet_ridge", "glmnet_elastic", "glmnet_lasso",
    "ranger_sqrt_node5", "ranger_mtry15_node5", "ranger_mtry25_node10",
    "lightgbm_d2_l7", "lightgbm_d3_l15", "lightgbm_d4_l31",
    "svm_radial_c1", "svm_radial_c4",
    "gbm_d2", "gbm_d3", "nnet4"
  ),
  top_n = c(NA_integer_, 20L, 10L),
  stringsAsFactors = FALSE
)

# Keep the first run broad but not excessive.
all_auc <- list()
all_pred <- list()
for (ii in seq_len(nrow(model_grid))) {
  tr <- model_grid$transform[ii]
  mn <- model_grid$method[ii]
  tn <- model_grid$top_n[ii]
  # SVM and neural network are expensive and unstable with too many redundant variants.
  if (mn %in% c("svm_radial_c1", "svm_radial_c4", "nnet4") && tr %in% c("raw", "rank")) next
  if (mn %in% c("lightgbm_d4_l31", "gbm_d3") && !is.na(tn) && tn == 10L) next
  x <- build_transform(cell_mat, tr)
  res <- run_cv_model(x, y, tr, mn, tn, n_repeats = 5, k = 5)
  all_auc[[length(all_auc) + 1L]] <- res$auc
  all_pred[[length(all_pred) + 1L]] <- res$predictions
}

auc_by_repeat <- do.call(rbind, all_auc)
predictions <- do.call(rbind, all_pred)

summary_df <- aggregate(cbind(AUC, Mean_R, Mean_NR) ~ model + transform + method + top_n,
                        data = auc_by_repeat, FUN = mean)
summary_df$MedianAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = median)$AUC[match(summary_df$model, aggregate(AUC ~ model, data = auc_by_repeat, FUN = median)$model)]
summary_df$SDAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = sd)$AUC[match(summary_df$model, aggregate(AUC ~ model, data = auc_by_repeat, FUN = sd)$model)]
summary_df$MaxAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = max)$AUC[match(summary_df$model, aggregate(AUC ~ model, data = auc_by_repeat, FUN = max)$model)]
summary_df$MinAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = min)$AUC[match(summary_df$model, aggregate(AUC ~ model, data = auc_by_repeat, FUN = min)$model)]
summary_df <- summary_df[order(summary_df$AUC, decreasing = TRUE), ]

fwrite(auc_by_repeat, file.path(out_dir, "cell_fraction_only_repeated5fold_auc_by_repeat.txt"), sep = "\t")
fwrite(predictions, file.path(out_dir, "cell_fraction_only_repeated5fold_predictions.txt"), sep = "\t")
fwrite(summary_df, file.path(out_dir, "cell_fraction_only_repeated5fold_auc_summary.txt"), sep = "\t")

# Exploratory OOF ensemble from the top cell-only models. This is a screening summary, not an external validation.
top_models <- head(summary_df$model, 10)
ens_rows <- list()
for (rr in sort(unique(predictions$repeat_id))) {
  wide <- predictions[predictions$repeat_id == rr & predictions$model %in% top_models, ]
  if (!nrow(wide)) next
  pred_mat <- tapply(wide$prediction, list(wide$sample, wide$model), mean)
  yy <- y[rownames(pred_mat)]
  for (nn in c(3, 5, 10)) {
    use_models <- head(top_models, nn)
    pp <- rowMeans(pred_mat[, use_models, drop = FALSE], na.rm = TRUE)
    ens_rows[[length(ens_rows) + 1L]] <- data.frame(
      model = paste0("OOF_mean_top", nn, "_cell_only_models"),
      repeat_id = rr,
      AUC = auc_manual(yy, pp),
      Mean_R = mean(pp[yy == 1], na.rm = TRUE),
      Mean_NR = mean(pp[yy == 0], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
ensemble_auc <- do.call(rbind, ens_rows)
fwrite(ensemble_auc, file.path(out_dir, "cell_fraction_only_exploratory_oof_ensemble_auc.txt"), sep = "\t")

best_model <- summary_df$model[1]
best_pred <- predictions[predictions$model == best_model & predictions$repeat_id == 1, ]
best_pred$response <- factor(ifelse(best_pred$y == 1, "Responder", "Non-responder"),
                             levels = c("Non-responder", "Responder"))

p1 <- ggplot(head(summary_df, 30), aes(x = reorder(model, AUC), y = AUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#F66463") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Mean repeated 5-fold CV AUC",
       title = "Cell-fraction-only model optimization") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cell_fraction_only_top30_auc_barplot.pdf"), p1, width = 9, height = 7)
ggsave(file.path(out_dir, "cell_fraction_only_top30_auc_barplot.png"), p1, width = 9, height = 7, dpi = 300)

p2 <- ggplot(best_pred, aes(x = response, y = prediction, fill = response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
  scale_fill_manual(values = c("Non-responder" = "#95A8AC", "Responder" = "#F66463")) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Predicted response probability",
       title = paste0("Best cell-only model: ", best_model),
       subtitle = paste0("Repeat 1 AUC = ", round(auc_manual(best_pred$y, best_pred$prediction), 3))) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cell_fraction_only_best_model_response_boxplot_repeat1.pdf"), p2, width = 4.2, height = 4.2)
ggsave(file.path(out_dir, "cell_fraction_only_best_model_response_boxplot_repeat1.png"), p2, width = 4.2, height = 4.2, dpi = 300)

writeLines(c(
  paste0("Best cell-fraction-only model: ", best_model),
  paste0("Mean AUC: ", signif(summary_df$AUC[1], 5)),
  paste0("Median AUC: ", signif(summary_df$MedianAUC[1], 5)),
  paste0("Max repeat AUC: ", signif(summary_df$MaxAUC[1], 5)),
  "",
  "Only ImmuCellAI2 cell fractions were used.",
  "No IFN-gamma, cytotoxicity, checkpoint, MHC, T cell inflamed GEP, cancer type, drug, target, study, or batch features were included.",
  "OOF ensemble results are exploratory and should not be treated as independent validation."
), con = file.path(out_dir, "cell_fraction_only_model_readme.txt"))

print(head(summary_df, 30))
print(ensemble_auc)
