options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
out_dir <- file.path(base_dir, "cell_fraction_only_model_optimization_fast")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(ranger)
  library(lightgbm)
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

make_stratified_folds <- function(y, k = 5L) {
  folds <- integer(length(y))
  for (yy in sort(unique(y))) {
    idx <- sample(which(y == yy))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

scale_by_train <- function(x_train, x_test) {
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
  if (transform_name == "asin") {
    z <- asin(sqrt(x))
  } else if (transform_name == "logit") {
    p <- pmin(pmax(x, eps), 1 - eps)
    z <- log(p / (1 - p))
  } else if (transform_name == "clr") {
    lx <- log(x + eps)
    z <- sweep(lx, 1, rowMeans(lx), "-")
  } else {
    stop("Unknown transform: ", transform_name)
  }
  colnames(z) <- paste0(transform_name, "_", make.names(colnames(x)))
  rownames(z) <- rownames(x)
  z
}

predict_glmnet <- function(y_train, x_train, x_test, alpha = 0) {
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(
      x_train, y_train, family = "binomial", alpha = alpha,
      type.measure = "auc", nfolds = 5
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

predict_ranger <- function(y_train, x_train, x_test, mtry = 15L, min_node = 5L) {
  dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train, check.names = FALSE))
  fit <- tryCatch(
    ranger(
      y ~ ., data = dat, probability = TRUE, classification = TRUE,
      num.trees = 800, mtry = min(mtry, ncol(x_train)),
      min.node.size = min_node, seed = sample.int(1e6, 1),
      num.threads = n_threads
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  pr <- predict(fit, data = as.data.frame(x_test, check.names = FALSE))$predictions
  as.numeric(pr[, "1"])
}

predict_lgb <- function(y_train, x_train, x_test, leaves = 15L, depth = 3L, lr = 0.03, rounds = 350L) {
  dtrain <- tryCatch(lgb.Dataset(data = as.matrix(x_train), label = y_train), error = function(e) NULL)
  if (is.null(dtrain)) return(rep(NA_real_, nrow(x_test)))
  params <- list(
    objective = "binary",
    metric = "auc",
    learning_rate = lr,
    num_leaves = leaves,
    max_depth = depth,
    min_data_in_leaf = 25,
    feature_fraction = 0.85,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    lambda_l1 = 0.05,
    lambda_l2 = 1,
    verbosity = -1,
    num_threads = n_threads,
    force_col_wise = TRUE
  )
  fit <- tryCatch(
    lgb.train(params = params, data = dtrain, nrounds = rounds, verbose = -1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, as.matrix(x_test)))
}

predict_model <- function(model_name, y_train, x_train, x_test) {
  if (model_name == "glmnet_ridge") return(predict_glmnet(y_train, x_train, x_test, alpha = 0))
  if (model_name == "glmnet_elastic") return(predict_glmnet(y_train, x_train, x_test, alpha = 0.5))
  if (model_name == "ranger_mtry15_node5") return(predict_ranger(y_train, x_train, x_test, mtry = 15, min_node = 5))
  if (model_name == "ranger_mtry25_node10") return(predict_ranger(y_train, x_train, x_test, mtry = 25, min_node = 10))
  if (model_name == "lightgbm_d2_l7") return(predict_lgb(y_train, x_train, x_test, leaves = 7, depth = 2, lr = 0.03, rounds = 350))
  if (model_name == "lightgbm_d3_l15") return(predict_lgb(y_train, x_train, x_test, leaves = 15, depth = 3, lr = 0.03, rounds = 350))
  stop("Unknown model: ", model_name)
}

run_one_model <- function(x_all, y, transform_name, model_name, top_n, n_repeats = 3L, k = 5L) {
  label <- paste(transform_name, model_name, ifelse(is.na(top_n), "all53", paste0("top", top_n)), sep = "__")
  auc_rows <- vector("list", n_repeats)
  pred_rows <- vector("list", n_repeats)
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
    auc_rows[[rr]] <- data.frame(
      model = label, transform = transform_name, method = model_name,
      top_n = ifelse(is.na(top_n), "all53", paste0("top", top_n)),
      repeat_id = rr, AUC = auc_manual(y, pred),
      Mean_R = mean(pred[y == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y == 0], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    pred_rows[[rr]] <- data.frame(
      model = label, repeat_id = rr, sample = names(y),
      y = y, prediction = pred, stringsAsFactors = FALSE
    )
    message(label, " repeat ", rr, "/", n_repeats, " AUC=", round(auc_rows[[rr]]$AUC, 4))
  }
  list(auc = do.call(rbind, auc_rows), pred = do.call(rbind, pred_rows))
}

message("Reading response and cell fractions...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
fraction_df <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
names(fraction_df)[1] <- "Run"
rownames(fraction_df) <- as.character(fraction_df$Run)
fraction_df$Run <- NULL

response_binary <- rep(NA_integer_, nrow(score_df))
response_binary[score_df$ResponseGroup == "Responder"] <- 1L
response_binary[score_df$ResponseGroup == "Non-responder"] <- 0L
response_binary[is.na(response_binary) & score_df$Response == "R"] <- 1L
response_binary[is.na(response_binary) & score_df$Response == "NR"] <- 0L
score_df$ResponseBinaryClean <- response_binary
score_df <- score_df[!is.na(score_df$ResponseBinaryClean), , drop = FALSE]
score_df$Run <- as.character(score_df$Run)
score_df <- score_df[score_df$Run %in% rownames(fraction_df), , drop = FALSE]
cell_mat <- as.matrix(fraction_df[score_df$Run, , drop = FALSE])
storage.mode(cell_mat) <- "numeric"
y <- as.integer(score_df$ResponseBinaryClean)
names(y) <- score_df$Run
message("Samples: ", length(y), "; R: ", sum(y == 1), "; NR: ", sum(y == 0), "; cells: ", ncol(cell_mat))

grid <- expand.grid(
  transform = c("logit", "asin", "clr"),
  method = c("glmnet_ridge", "glmnet_elastic", "ranger_mtry15_node5",
             "ranger_mtry25_node10", "lightgbm_d2_l7", "lightgbm_d3_l15"),
  top_n = c(NA_integer_, 20L),
  stringsAsFactors = FALSE
)

auc_file <- file.path(out_dir, "cell_fraction_only_fast_repeated5fold_auc_by_repeat.txt")
pred_file <- file.path(out_dir, "cell_fraction_only_fast_repeated5fold_predictions.txt")
if (file.exists(auc_file)) file.remove(auc_file)
if (file.exists(pred_file)) file.remove(pred_file)

all_auc <- list()
for (ii in seq_len(nrow(grid))) {
  tr <- grid$transform[ii]
  mn <- grid$method[ii]
  tn <- grid$top_n[ii]
  x <- build_transform(cell_mat, tr)
  res <- run_one_model(x, y, tr, mn, tn, n_repeats = 3, k = 5)
  all_auc[[length(all_auc) + 1L]] <- res$auc
  fwrite(res$auc, auc_file, sep = "\t", append = file.exists(auc_file))
  fwrite(res$pred, pred_file, sep = "\t", append = file.exists(pred_file))
}

auc_by_repeat <- do.call(rbind, all_auc)
summary_df <- aggregate(cbind(AUC, Mean_R, Mean_NR) ~ model + transform + method + top_n,
                        data = auc_by_repeat, FUN = mean)
med <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = median)
sdv <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = sd)
mx <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = max)
mn <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = min)
summary_df$MedianAUC <- med$AUC[match(summary_df$model, med$model)]
summary_df$SDAUC <- sdv$AUC[match(summary_df$model, sdv$model)]
summary_df$MaxAUC <- mx$AUC[match(summary_df$model, mx$model)]
summary_df$MinAUC <- mn$AUC[match(summary_df$model, mn$model)]
summary_df <- summary_df[order(summary_df$AUC, decreasing = TRUE), ]
fwrite(summary_df, file.path(out_dir, "cell_fraction_only_fast_auc_summary.txt"), sep = "\t")

p <- ggplot(head(summary_df, 25), aes(x = reorder(model, AUC), y = AUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Mean repeated 5-fold CV AUC",
       title = "Cell-fraction-only model screening") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "cell_fraction_only_fast_top25_auc_barplot.pdf"), p, width = 9, height = 6)
ggsave(file.path(out_dir, "cell_fraction_only_fast_top25_auc_barplot.png"), p, width = 9, height = 6, dpi = 300)

writeLines(c(
  paste0("Best model: ", summary_df$model[1]),
  paste0("Mean AUC: ", signif(summary_df$AUC[1], 5)),
  paste0("Median AUC: ", signif(summary_df$MedianAUC[1], 5)),
  paste0("Max repeat AUC: ", signif(summary_df$MaxAUC[1], 5)),
  "Only ImmuCellAI2 cell fraction features were used.",
  "No gene-expression signatures or clinical context variables were included."
), file.path(out_dir, "cell_fraction_only_fast_readme.txt"))

print(head(summary_df, 25))
