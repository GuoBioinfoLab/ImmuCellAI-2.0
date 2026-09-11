options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
axis_dir <- file.path(base_dir, "C3_C4_axis_ICB_response_scores")
out_dir <- file.path(base_dir, "unified_elasticnet_ICB_score_LOSO")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
axis_score_file <- file.path(axis_dir, "external_ICB_C3_C4_axis_scores.txt")
meta_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
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

row_mean_existing <- function(mat, cells) {
  cells <- intersect(cells, colnames(mat))
  if (!length(cells)) return(rep(NA_real_, nrow(mat)))
  rowMeans(mat[, cells, drop = FALSE], na.rm = TRUE)
}

make_feature_name_safe <- function(x) make.names(x, unique = TRUE)

scale_train_test <- function(x_train, x_test) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  x_train[!is.finite(x_train)] <- NA_real_
  x_test[!is.finite(x_test)] <- NA_real_
  center <- apply(x_train, 2, median, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  for (j in seq_len(ncol(x_train))) {
    x_train[!is.finite(x_train[, j]), j] <- center[j]
    x_test[!is.finite(x_test[, j]), j] <- center[j]
  }
  scalev <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(scalev) & scalev > 0
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

fit_predict_glmnet <- function(x_train, y_train, x_test, alpha) {
  weights <- ifelse(y_train == 1, 0.5 / mean(y_train == 1), 0.5 / mean(y_train == 0))
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(
      x_train, y_train,
      family = "binomial",
      alpha = alpha,
      weights = weights,
      type.measure = "auc",
      nfolds = 5,
      standardize = FALSE
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(pred = rep(NA_real_, nrow(x_test)), fit = NULL))
  pred <- as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
  list(pred = pred, fit = fit)
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

message("Samples: ", length(y), "; responders: ", sum(y == 1), "; nonresponders: ", sum(y == 0))

raw <- fraction
asin <- asin(sqrt(pmax(raw, 0)))
clr <- clr_transform(raw)
logit <- logit_transform(raw)
colnames(raw) <- paste0("raw_", make_feature_name_safe(colnames(raw)))
colnames(asin) <- paste0("asin_", make_feature_name_safe(gsub("^raw_", "", colnames(raw))))
colnames(clr) <- paste0("clr_", make_feature_name_safe(gsub("^raw_", "", colnames(raw))))
colnames(logit) <- paste0("logit_", make_feature_name_safe(gsub("^raw_", "", colnames(raw))))

axis_keep <- c(
  "C1_spearman_raw", "C2_spearman_raw", "C3_spearman_raw", "C4_spearman_raw",
  "C1_pearson_clr", "C2_pearson_clr", "C3_pearson_clr", "C4_pearson_clr",
  "C4_minus_C3_spearman_raw", "C4_minus_C3_pearson_clr", "C4_minus_C3_cosine_raw",
  "auto_top5_C4sum_over_C3sum_logratio", "auto_top10_C4sum_over_C3sum_logratio",
  "manual_effector_minus_C3suppressive_raw",
  "manual_effector_over_C3suppressive_logratio",
  "C3_MDSC_M1_negative",
  "C4effector_over_MDSC_M1_logratio",
  "TCGA_C4_vs_C3_glmnet_ridge_probability_asin",
  "TCGA_C4_vs_C3_glmnet_lasso_probability_asin",
  "TCGA_C4_vs_C3_ranger_probability_asin",
  "TCGA_C4_vs_C3_ranger_probability_clr"
)
axis_keep <- intersect(axis_keep, colnames(axis_scores))
axis_mat <- as.matrix(axis_scores[, axis_keep, drop = FALSE])
storage.mode(axis_mat) <- "numeric"
colnames(axis_mat) <- paste0("axis_", make_feature_name_safe(colnames(axis_mat)))

cell_names <- gsub("^raw_", "", colnames(raw))
raw_for_bio <- fraction
effector_cells <- c("CD8Tem", "CD8Temra", "CD8Trm", "Tc", "Tex", "CD8Tcm", "CD8Tn",
                    "cNK", "NKreg", "NKT", "gdT", "MAIT", "Th1", "Tfh",
                    "CD4Tem", "CD4Temra", "CD4Trm", "cDC1", "cDC2", "pDC", "PB", "PC", "MBC", "BGC")
suppressive_cells <- c("MDSC", "M1", "M2", "TAM", "cMo", "intMo", "ncMo", "Neutrophil", "Treg")
myeloid_cells <- c("MDSC", "M1", "M2", "TAM", "cMo", "intMo", "ncMo", "Neutrophil", "M0")
tcell_cells <- c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm", "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex", "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg")

bio <- data.frame(
  effector_sum = row_sum_existing(raw_for_bio, effector_cells),
  suppressive_sum = row_sum_existing(raw_for_bio, suppressive_cells),
  myeloid_sum = row_sum_existing(raw_for_bio, myeloid_cells),
  tcell_sum = row_sum_existing(raw_for_bio, tcell_cells),
  MDSC = row_sum_existing(raw_for_bio, "MDSC"),
  M1 = row_sum_existing(raw_for_bio, "M1"),
  TAM = row_sum_existing(raw_for_bio, "TAM"),
  Treg = row_sum_existing(raw_for_bio, "Treg"),
  effector_over_suppressive = log((row_sum_existing(raw_for_bio, effector_cells) + 1e-5) /
                                    (row_sum_existing(raw_for_bio, suppressive_cells) + 1e-5)),
  tcell_over_myeloid = log((row_sum_existing(raw_for_bio, tcell_cells) + 1e-5) /
                             (row_sum_existing(raw_for_bio, myeloid_cells) + 1e-5)),
  stringsAsFactors = FALSE
)
colnames(bio) <- paste0("bio_", make_feature_name_safe(colnames(bio)))
bio <- as.matrix(bio)

feature_sets <- list(
  axis_only = axis_mat,
  biology_compact = bio,
  cell_asin53 = asin,
  cell_clr53 = clr,
  cell_logit53 = logit,
  axis_plus_biology = cbind(axis_mat, bio),
  axis_plus_clr53 = cbind(axis_mat, clr),
  axis_plus_asin53 = cbind(axis_mat, asin),
  axis_bio_clr53 = cbind(axis_mat, bio, clr),
  all_cell_transforms_axis_bio = cbind(axis_mat, bio, asin, clr, logit)
)

valid_studies <- names(which(table(meta$SRA_study) >= 20))
valid_studies <- valid_studies[sapply(valid_studies, function(st) {
  yy <- y[meta$SRA_study == st]
  sum(yy == 1) >= 5 && sum(yy == 0) >= 5
})]
message("Valid LOSO studies: ", length(valid_studies))

alphas <- c(0, 0.25, 0.5, 0.75, 1)
summary_rows <- list()
pred_rows <- list()

for (fs_name in names(feature_sets)) {
  x_all <- as.matrix(feature_sets[[fs_name]])
  for (alpha in alphas) {
    model_name <- paste0(fs_name, "__elasticnet_alpha", alpha)
    message("Running ", model_name)
    study_rows <- list()
    for (st in valid_studies) {
      test_idx <- which(meta$SRA_study == st)
      train_idx <- which(meta$SRA_study != st)
      xs <- scale_train_test(x_all[train_idx, , drop = FALSE], x_all[test_idx, , drop = FALSE])
      pred <- fit_predict_glmnet(xs$train, y[train_idx], xs$test, alpha = alpha)$pred
      a <- auc_manual(y[test_idx], pred)
      study_rows[[length(study_rows) + 1L]] <- data.frame(
        model = model_name,
        feature_set = fs_name,
        alpha = alpha,
        SRA_study = st,
        N = length(test_idx),
        N_R = sum(y[test_idx] == 1),
        N_NR = sum(y[test_idx] == 0),
        AUC = a,
        Mean_R = mean(pred[y[test_idx] == 1], na.rm = TRUE),
        Mean_NR = mean(pred[y[test_idx] == 0], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        model = model_name,
        sample = names(y)[test_idx],
        SRA_study = st,
        y = y[test_idx],
        prediction = pred,
        stringsAsFactors = FALSE
      )
    }
    st_df <- rbindlist(study_rows)
    summary_rows[[length(summary_rows) + 1L]] <- st_df[, .(
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
    ), by = .(model, feature_set, alpha)]
    fwrite(st_df, file.path(out_dir, "LOSO_elasticnet_auc_by_study_running.txt"),
           sep = "\t", append = file.exists(file.path(out_dir, "LOSO_elasticnet_auc_by_study_running.txt")))
  }
}

summary_df <- rbindlist(summary_rows)
summary_df <- summary_df[order(-N_AUC_ge_0.75, -MeanAUC, -WeightedAUC)]
pred_df <- rbindlist(pred_rows)
auc_by_study <- fread(file.path(out_dir, "LOSO_elasticnet_auc_by_study_running.txt"), data.table = FALSE)

fwrite(summary_df, file.path(out_dir, "LOSO_elasticnet_model_summary.txt"), sep = "\t")
fwrite(auc_by_study, file.path(out_dir, "LOSO_elasticnet_auc_by_study.txt"), sep = "\t")
fwrite(pred_df, file.path(out_dir, "LOSO_elasticnet_predictions.txt"), sep = "\t")

best_model <- summary_df$model[1]
best_auc <- auc_by_study[auc_by_study$model == best_model, ]
best_pred <- pred_df[pred_df$model == best_model, ]
fwrite(best_auc, file.path(out_dir, "BEST_LOSO_elasticnet_auc_by_study.txt"), sep = "\t")
fwrite(best_pred, file.path(out_dir, "BEST_LOSO_elasticnet_predictions.txt"), sep = "\t")

# Fit final unified score on all ICB samples for downstream scoring and coefficient interpretation.
best_feature_set <- summary_df$feature_set[1]
best_alpha <- summary_df$alpha[1]
x_final <- as.matrix(feature_sets[[best_feature_set]])
xs_final <- scale_train_test(x_final, x_final)
final_fit <- fit_predict_glmnet(xs_final$train, y, xs_final$test, alpha = best_alpha)$fit
final_score <- as.numeric(predict(final_fit, newx = xs_final$test, s = "lambda.min", type = "response"))
coef_mat <- as.matrix(coef(final_fit, s = "lambda.min"))
coef_df <- data.frame(feature = rownames(coef_mat), coefficient = as.numeric(coef_mat[, 1]), row.names = NULL)
coef_df <- coef_df[coef_df$coefficient != 0, ]
coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]

final_score_df <- data.frame(
  Run = meta$Run,
  SRA_study = meta$SRA_study,
  Cancer = meta$Cancer,
  Cancer_type = meta$Cancer_type,
  Drug = meta$Drug,
  disease = meta$disease,
  Response = meta$Response,
  response_binary = y,
  UnifiedElasticNetICBScore = final_score,
  stringsAsFactors = FALSE
)
fwrite(final_score_df, file.path(out_dir, "final_unified_elasticnet_ICB_score_trained_all_samples.txt"), sep = "\t")
fwrite(coef_df, file.path(out_dir, "final_unified_elasticnet_nonzero_coefficients.txt"), sep = "\t")

p1 <- ggplot(head(summary_df, 25), aes(x = reorder(model, MeanAUC), y = MeanAUC, fill = N_AUC_ge_0.75)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  coord_flip() +
  scale_fill_gradient(low = "#95A8AC", high = "#F66463") +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = "Mean LOSO AUC", fill = "N AUC >= 0.75",
       title = "Unified elastic-net ICB score model screening")
ggsave(file.path(out_dir, "LOSO_elasticnet_model_summary_top25.pdf"), p1, width = 9, height = 6)
ggsave(file.path(out_dir, "LOSO_elasticnet_model_summary_top25.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(best_auc, aes(x = reorder(SRA_study, AUC), y = AUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.75, color = "#F66463", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "grey40", linetype = "dotted") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Leave-one-study-out AUC",
       title = paste0("Best unified elastic-net ICB score: ", best_model))
ggsave(file.path(out_dir, "BEST_LOSO_elasticnet_auc_by_study.pdf"), p2, width = 7.5, height = 5.5)
ggsave(file.path(out_dir, "BEST_LOSO_elasticnet_auc_by_study.png"), p2, width = 7.5, height = 5.5, dpi = 300)

writeLines(c(
  paste0("Best model: ", best_model),
  paste0("Feature set: ", best_feature_set),
  paste0("Alpha: ", best_alpha),
  paste0("Mean LOSO AUC: ", signif(summary_df$MeanAUC[1], 5)),
  paste0("Weighted LOSO AUC: ", signif(summary_df$WeightedAUC[1], 5)),
  paste0("Number of studies with AUC >= 0.75: ", summary_df$N_AUC_ge_0.75[1], "/", summary_df$N_studies[1]),
  paste0("Number of studies with AUC >= 0.80: ", summary_df$N_AUC_ge_0.80[1], "/", summary_df$N_studies[1]),
  "Predictors use only ImmuCellAI2 cell fractions and C3/C4 axis features derived from cell fractions.",
  "No IFN-gamma, cytotoxicity, checkpoint, MHC, T cell-inflamed GEP, cancer type, drug, target, study, or batch predictors were used."
), file.path(out_dir, "README_unified_elasticnet_ICB_score_LOSO.txt"))

print(head(summary_df, 20))
print(best_auc[order(-best_auc$AUC), ])
print(head(coef_df, 30))
