options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
axis_dir <- file.path(base_dir, "C3_C4_axis_ICB_response_scores")
more_model_dir <- file.path(base_dir, "more_unified_ICB_score_models")
out_dir <- file.path(base_dir, "score_cluster_joint_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
axis_score_file <- file.path(axis_dir, "external_ICB_C3_C4_axis_scores.txt")
meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
previous_loso_prediction_file <- file.path(more_model_dir, "LOSO_predictions.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ranger)
  library(lightgbm)
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

predict_ranger_prob <- function(x_train, y_train, x_test, mtry = NULL, min_node = 8) {
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(ncol(x_train))))
  cw <- class_weights(y_train)
  fit <- ranger(
    x = as.data.frame(x_train, check.names = FALSE),
    y = factor(y_train, levels = c(0, 1)),
    probability = TRUE,
    classification = TRUE,
    num.trees = 1200,
    mtry = min(mtry, ncol(x_train)),
    min.node.size = min_node,
    class.weights = c("0" = cw[which(y_train == 0)[1]], "1" = cw[which(y_train == 1)[1]]),
    seed = sample.int(1e6, 1),
    num.threads = n_threads
  )
  list(
    train = as.numeric(predict(fit, data = as.data.frame(x_train, check.names = FALSE))$predictions[, "1"]),
    test = as.numeric(predict(fit, data = as.data.frame(x_test, check.names = FALSE))$predictions[, "1"])
  )
}

predict_glmnet_prob <- function(x_train, y_train, x_test, alpha = 0) {
  fit <- suppressWarnings(cv.glmnet(
    x_train, y_train,
    family = "binomial",
    alpha = alpha,
    weights = class_weights(y_train),
    type.measure = "auc",
    nfolds = 5,
    standardize = FALSE
  ))
  as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

predict_lightgbm_prob <- function(x_train, y_train, x_test) {
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
  as.numeric(predict(fit, x_test))
}

run_loso_model <- function(x, y, study, valid_studies, model_name) {
  rows <- list()
  preds <- list()
  for (st in valid_studies) {
    test_idx <- which(study == st)
    train_idx <- which(study != st)
    xs <- scale_train_test(x[train_idx, , drop = FALSE], x[test_idx, , drop = FALSE])
    pred <- tryCatch({
      if (model_name == "glmnet_ridge") {
        predict_glmnet_prob(xs$train, y[train_idx], xs$test, alpha = 0)
      } else if (model_name == "glmnet_elastic") {
        predict_glmnet_prob(xs$train, y[train_idx], xs$test, alpha = 0.5)
      } else if (model_name == "ranger_balanced") {
        predict_ranger_prob(xs$train, y[train_idx], xs$test)$test
      } else if (model_name == "lightgbm_shallow") {
        predict_lightgbm_prob(xs$train, y[train_idx], xs$test)
      } else {
        stop("Unknown model")
      }
    }, error = function(e) rep(NA_real_, length(test_idx)))

    rows[[length(rows) + 1L]] <- data.frame(
      SRA_study = st,
      N = length(test_idx),
      N_R = sum(y[test_idx] == 1),
      N_NR = sum(y[test_idx] == 0),
      AUC = auc_manual(y[test_idx], pred),
      Mean_R = mean(pred[y[test_idx] == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y[test_idx] == 0], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    preds[[length(preds) + 1L]] <- data.frame(
      sample = names(y)[test_idx],
      SRA_study = st,
      y = y[test_idx],
      prediction = pred,
      stringsAsFactors = FALSE
    )
  }
  list(auc = rbindlist(rows), pred = rbindlist(preds))
}

run_stacked_score_cluster_loso <- function(clr_x, cluster_x, y, study, valid_studies) {
  rows <- list()
  preds <- list()
  for (st in valid_studies) {
    test_idx <- which(study == st)
    train_idx <- which(study != st)
    clr_scaled <- scale_train_test(clr_x[train_idx, , drop = FALSE], clr_x[test_idx, , drop = FALSE])
    base <- predict_ranger_prob(clr_scaled$train, y[train_idx], clr_scaled$test)

    meta_train <- cbind(RF_CLR_score = base$train, cluster_x[train_idx, , drop = FALSE])
    meta_test <- cbind(RF_CLR_score = base$test, cluster_x[test_idx, , drop = FALSE])
    meta_scaled <- scale_train_test(meta_train, meta_test)
    pred <- tryCatch(predict_glmnet_prob(meta_scaled$train, y[train_idx], meta_scaled$test, alpha = 0),
                     error = function(e) rep(NA_real_, length(test_idx)))

    rows[[length(rows) + 1L]] <- data.frame(
      SRA_study = st,
      N = length(test_idx),
      N_R = sum(y[test_idx] == 1),
      N_NR = sum(y[test_idx] == 0),
      AUC = auc_manual(y[test_idx], pred),
      Mean_R = mean(pred[y[test_idx] == 1], na.rm = TRUE),
      Mean_NR = mean(pred[y[test_idx] == 0], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    preds[[length(preds) + 1L]] <- data.frame(
      sample = names(y)[test_idx],
      SRA_study = st,
      y = y[test_idx],
      prediction = pred,
      stringsAsFactors = FALSE
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

clr <- clr_transform(fraction)
colnames(clr) <- paste0("clr_", make.names(colnames(fraction), unique = TRUE))

c_spear <- as.matrix(axis_scores[, paste0("C", 1:4, "_spearman_raw"), drop = FALSE])
colnames(c_spear) <- paste0("cluster_", colnames(c_spear))
nearest_cluster_spearman <- paste0("C", apply(c_spear, 1, which.max))

c_clr <- as.matrix(axis_scores[, paste0("C", 1:4, "_pearson_clr"), drop = FALSE])
colnames(c_clr) <- paste0("cluster_", colnames(c_clr))
nearest_cluster_clr <- paste0("C", apply(c_clr, 1, which.max))

cluster_onehot <- model.matrix(~ nearest_cluster_spearman - 1)
colnames(cluster_onehot) <- gsub("^nearest_cluster_spearman", "nearest_", colnames(cluster_onehot))

cluster_axis_features <- data.frame(
  C1_spearman = axis_scores$C1_spearman_raw,
  C2_spearman = axis_scores$C2_spearman_raw,
  C3_spearman = axis_scores$C3_spearman_raw,
  C4_spearman = axis_scores$C4_spearman_raw,
  C1_pearson_clr = axis_scores$C1_pearson_clr,
  C2_pearson_clr = axis_scores$C2_pearson_clr,
  C3_pearson_clr = axis_scores$C3_pearson_clr,
  C4_pearson_clr = axis_scores$C4_pearson_clr,
  C3_low_clr = -axis_scores$C3_pearson_clr,
  C4_minus_C3_spearman = axis_scores$C4_minus_C3_spearman_raw,
  C4_minus_C3_pearson_clr = axis_scores$C4_minus_C3_pearson_clr,
  C4_minus_C3_cosine = axis_scores$C4_minus_C3_cosine_raw,
  stringsAsFactors = FALSE
)
cluster_axis_features <- as.matrix(cluster_axis_features)
colnames(cluster_axis_features) <- paste0("axis_", make.names(colnames(cluster_axis_features), unique = TRUE))

cluster_features <- cbind(cluster_axis_features, cluster_onehot)
clr_plus_cluster <- cbind(clr, cluster_features)

valid_studies <- names(which(table(study) >= 20))
valid_studies <- valid_studies[sapply(valid_studies, function(st) {
  yy <- y[study == st]
  sum(yy == 1) >= 5 && sum(yy == 0) >= 5
})]

message("Valid studies: ", length(valid_studies))

# Descriptive link between RF score and cluster states.
previous_pred <- fread(previous_loso_prediction_file, data.table = FALSE)
rf_pred <- previous_pred[previous_pred$model == "cell_clr53__ranger_balanced", , drop = FALSE]
rf_pred <- rf_pred[match(meta$Run, rf_pred$sample), , drop = FALSE]

cluster_map <- data.frame(
  Run = meta$Run,
  SRA_study = study,
  Response = meta$Response,
  response_binary = y,
  nearest_cluster_spearman = nearest_cluster_spearman,
  nearest_cluster_clr = nearest_cluster_clr,
  RF_CLR_LOSO_score = rf_pred$prediction,
  C3_low_clr = -axis_scores$C3_pearson_clr,
  C4_minus_C3_spearman = axis_scores$C4_minus_C3_spearman_raw,
  stringsAsFactors = FALSE
)
fwrite(cluster_map, file.path(out_dir, "ICB_RF_score_with_nearest_cluster.txt"), sep = "\t")

cluster_summary <- as.data.table(cluster_map)[, .(
  N = .N,
  N_R = sum(response_binary == 1),
  ResponseRate = mean(response_binary == 1),
  Mean_RF_CLR_LOSO_score = mean(RF_CLR_LOSO_score, na.rm = TRUE),
  Median_RF_CLR_LOSO_score = median(RF_CLR_LOSO_score, na.rm = TRUE),
  AUC_RF_within_cluster = auc_manual(response_binary, RF_CLR_LOSO_score),
  Mean_C3_low_clr = mean(C3_low_clr, na.rm = TRUE)
), by = nearest_cluster_spearman][order(nearest_cluster_spearman)]
fwrite(cluster_summary, file.path(out_dir, "cluster_summary_RF_score_response.txt"), sep = "\t")

model_results <- list()
model_predictions <- list()

configs <- list(
  cluster_axis_glmnet_ridge = list(x = cluster_features, model = "glmnet_ridge"),
  cluster_axis_ranger = list(x = cluster_features, model = "ranger_balanced"),
  clr_plus_cluster_glmnet_ridge = list(x = clr_plus_cluster, model = "glmnet_ridge"),
  clr_plus_cluster_ranger = list(x = clr_plus_cluster, model = "ranger_balanced"),
  clr_plus_cluster_lightgbm = list(x = clr_plus_cluster, model = "lightgbm_shallow")
)

for (nm in names(configs)) {
  message("Running ", nm)
  res <- run_loso_model(configs[[nm]]$x, y, study, valid_studies, configs[[nm]]$model)
  res$auc$model <- nm
  res$pred$model <- nm
  model_results[[nm]] <- res$auc
  model_predictions[[nm]] <- res$pred
}

message("Running stacked_RF_CLR_score_plus_cluster_glmnet")
stacked <- run_stacked_score_cluster_loso(clr, cluster_features, y, study, valid_studies)
stacked$auc$model <- "stacked_RF_CLR_score_plus_cluster_glmnet"
stacked$pred$model <- "stacked_RF_CLR_score_plus_cluster_glmnet"
model_results[[stacked$auc$model[1]]] <- stacked$auc
model_predictions[[stacked$pred$model[1]]] <- stacked$pred

auc_by_study <- rbindlist(model_results, fill = TRUE)
pred_all <- rbindlist(model_predictions, fill = TRUE)

summary <- auc_by_study[, .(
  MeanAUC = mean(AUC, na.rm = TRUE),
  MedianAUC = median(AUC, na.rm = TRUE),
  WeightedAUC = weighted.mean(AUC, N, na.rm = TRUE),
  MinAUC = min(AUC, na.rm = TRUE),
  MaxAUC = max(AUC, na.rm = TRUE),
  N_studies = .N,
  N_AUC_ge_0.65 = sum(AUC >= 0.65, na.rm = TRUE),
  N_AUC_ge_0.70 = sum(AUC >= 0.70, na.rm = TRUE),
  N_AUC_ge_0.75 = sum(AUC >= 0.75, na.rm = TRUE),
  N_AUC_ge_0.80 = sum(AUC >= 0.80, na.rm = TRUE)
), by = model][order(-N_AUC_ge_0.75, -MeanAUC, -WeightedAUC)]

fwrite(auc_by_study, file.path(out_dir, "joint_score_cluster_LOSO_auc_by_study.txt"), sep = "\t")
fwrite(pred_all, file.path(out_dir, "joint_score_cluster_LOSO_predictions.txt"), sep = "\t")
fwrite(summary, file.path(out_dir, "joint_score_cluster_model_summary.txt"), sep = "\t")

p1 <- ggplot(cluster_map, aes(x = nearest_cluster_spearman, y = RF_CLR_LOSO_score, fill = factor(response_binary))) +
  geom_boxplot(outlier.shape = NA, color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.18, size = 0.6, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "#95A8AC", "1" = "#F66463"), labels = c("NR", "R")) +
  theme_bw(base_size = 11) +
  labs(x = "Nearest TCGA immune cluster", y = "RF-CLR LOSO ICB score", fill = "Response",
       title = "RF-CLR score distribution across C1-C4-like states") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "RF_CLR_score_by_nearest_cluster.pdf"), p1, width = 5.6, height = 4)
ggsave(file.path(out_dir, "RF_CLR_score_by_nearest_cluster.png"), p1, width = 5.6, height = 4, dpi = 300)

p2 <- ggplot(summary, aes(x = reorder(model, MeanAUC), y = MeanAUC, fill = N_AUC_ge_0.75)) +
  geom_col(color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  scale_fill_gradient(low = "#95A8AC", high = "#F66463") +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Mean LOSO AUC", fill = "N AUC >= 0.75",
       title = "Joint analysis of RF score and C1-C4 cluster features") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "joint_score_cluster_model_summary.pdf"), p2, width = 7.2, height = 4.8)
ggsave(file.path(out_dir, "joint_score_cluster_model_summary.png"), p2, width = 7.2, height = 4.8, dpi = 300)

best_model <- summary$model[1]
p3_df <- auc_by_study[auc_by_study$model == best_model, ][order(-AUC)]
p3 <- ggplot(p3_df, aes(x = reorder(SRA_study, AUC), y = AUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "grey40", linetype = "dotted") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "LOSO AUC", title = paste0("Best joint score-cluster model: ", best_model)) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "best_joint_score_cluster_auc_by_study.pdf"), p3, width = 7.2, height = 5.2)
ggsave(file.path(out_dir, "best_joint_score_cluster_auc_by_study.png"), p3, width = 7.2, height = 5.2, dpi = 300)

writeLines(c(
  "Joint analysis of RF-CLR ICB score and C1-C4 immune clusters",
  paste0("Best joint model: ", best_model),
  paste0("Mean LOSO AUC: ", signif(summary$MeanAUC[1], 5)),
  paste0("Weighted LOSO AUC: ", signif(summary$WeightedAUC[1], 5)),
  paste0("Number of studies with AUC >= 0.75: ", summary$N_AUC_ge_0.75[1], "/", summary$N_studies[1]),
  paste0("Number of studies with AUC >= 0.80: ", summary$N_AUC_ge_0.80[1], "/", summary$N_studies[1]),
  "Cluster labels were assigned by maximum similarity to TCGA C1-C4 centroids.",
  "No cancer type, drug, treatment target, study, batch, or gene-expression signature predictors were used."
), file.path(out_dir, "README_score_cluster_joint_analysis.txt"))

print(cluster_summary)
print(summary)
print(p3_df)
