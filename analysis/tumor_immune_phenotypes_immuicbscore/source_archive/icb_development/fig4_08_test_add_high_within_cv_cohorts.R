suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
out_dir <- file.path(base_dir, "fivefold_cv_unified_icb_score")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fraction_file <- file.path(
  fig4_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
clinical_file <- file.path(
  fig4_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "external_immunotherapy_C3_C4_scores_with_clinical.txt"
)
eligible_study_file <- file.path(
  fig4_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "more_unified_ICB_score_models",
  "BEST_by_LOSO_auc_by_study.txt"
)

original_high_loso <- c(
  "anti-PD1_SRP070710",
  "anti-PD1_SRP230414",
  "anti-PD1_SRP351936",
  "anti-PD1_ERP105482",
  "anti-PD1-anti-CTLA4_ERP105482"
)

added_high_within <- c(
  "anti-PD1_ERP107734",
  "anti-PD1_ERP105482",
  "anti-PD1_ERP117672"
)

expanded_studies <- unique(c(original_high_loso, added_high_within))

clr_transform <- function(x, eps = 1e-5) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx), "-")
}

fit_scaler <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  scalev <- apply(x, 2, sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  list(center = center, scale = scalev)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

auc_manual <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2) return(NA_real_)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

make_stratified_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  folds <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

train_predict_cv <- function(dat, cell_cols, folds, seed = 123, num_trees = 1200) {
  pred_rows <- vector("list", max(folds))
  for (fold in sort(unique(folds))) {
    train_dat <- dat[folds != fold, , drop = FALSE]
    valid_dat <- dat[folds == fold, , drop = FALSE]

    x_train_clr <- clr_transform(train_dat[, cell_cols, drop = FALSE])
    x_valid_clr <- clr_transform(valid_dat[, cell_cols, drop = FALSE])
    scaler <- fit_scaler(x_train_clr)
    x_train <- apply_scaler(x_train_clr, scaler)
    x_valid <- apply_scaler(x_valid_clr, scaler)

    y_train <- train_dat$ResponseBinary
    class_weights <- c(
      "0" = 0.5 / mean(y_train == 0L),
      "1" = 0.5 / mean(y_train == 1L)
    )

    set.seed(seed + fold)
    fit <- ranger(
      x = as.data.frame(x_train, check.names = FALSE),
      y = factor(y_train, levels = c(0, 1)),
      probability = TRUE,
      classification = TRUE,
      num.trees = num_trees,
      mtry = max(1, floor(sqrt(ncol(x_train)))),
      min.node.size = 8,
      class.weights = class_weights,
      seed = seed + fold,
      num.threads = 8
    )
    pred <- as.numeric(predict(fit, data = as.data.frame(x_valid, check.names = FALSE))$predictions[, "1"])
    pred_rows[[fold]] <- data.frame(
      sample = valid_dat$sample,
      SRA_study = valid_dat$SRA_study,
      ResponseBinary = valid_dat$ResponseBinary,
      fold = fold,
      Unified_ICB_Score = pred,
      stringsAsFactors = FALSE
    )
  }
  rbindlist(pred_rows, fill = TRUE)
}

summarize_pred <- function(pred, label, studies) {
  data.table(
    CV_design = label,
    Studies = paste(studies, collapse = ";"),
    N_studies = length(studies),
    N = nrow(pred),
    N_R = sum(pred$ResponseBinary == 1),
    N_NR = sum(pred$ResponseBinary == 0),
    AUC = auc_manual(pred$ResponseBinary, pred$Unified_ICB_Score),
    Mean_R = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 1], na.rm = TRUE),
    Mean_NR = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 0], na.rm = TRUE)
  )
}

fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
colnames(fraction)[1] <- "sample"
clinical <- fread(clinical_file, data.table = FALSE, check.names = TRUE)
eligible_auc <- fread(eligible_study_file, data.table = FALSE)

clinical_keep <- clinical[, c("Run", "SRA_study", "ResponseBinary", "ResponseGroup"), drop = FALSE]
colnames(clinical_keep)[1] <- "sample"
dat <- merge(fraction, clinical_keep, by = "sample")
dat <- dat[!is.na(dat$ResponseBinary) & dat$ResponseBinary %in% c(0, 1), ]
dat <- dat[!is.na(dat$SRA_study) & dat$SRA_study %in% unique(eligible_auc$SRA_study), ]
cell_cols <- setdiff(colnames(fraction), "sample")
dat <- dat[complete.cases(dat[, cell_cols, drop = FALSE]), ]
dat$ResponseBinary <- as.integer(dat$ResponseBinary)

dat_original <- dat[dat$SRA_study %in% original_high_loso, , drop = FALSE]
dat_expanded <- dat[dat$SRA_study %in% expanded_studies, , drop = FALSE]

folds_original <- make_stratified_folds(dat_original$ResponseBinary, k = 5, seed = 20260629)
folds_expanded <- make_stratified_folds(dat_expanded$ResponseBinary, k = 5, seed = 20260630)

pred_original <- train_predict_cv(dat_original, cell_cols, folds_original, seed = 7000)
pred_expanded <- train_predict_cv(dat_expanded, cell_cols, folds_expanded, seed = 8000)

summary_dt <- rbindlist(list(
  summarize_pred(pred_original, "Original high-LOSO 5 studies: pooled 5-fold CV", original_high_loso),
  summarize_pred(pred_expanded, "Expanded 7 studies: add high within-study CV cohorts", expanded_studies)
), fill = TRUE)

per_study_expanded <- as.data.table(pred_expanded)[, .(
  N = .N,
  N_R = sum(ResponseBinary == 1),
  N_NR = sum(ResponseBinary == 0),
  AUC = auc_manual(ResponseBinary, Unified_ICB_Score),
  Mean_R = mean(Unified_ICB_Score[ResponseBinary == 1], na.rm = TRUE),
  Mean_NR = mean(Unified_ICB_Score[ResponseBinary == 0], na.rm = TRUE)
), by = SRA_study][order(-AUC)]

fwrite(summary_dt, file.path(out_dir, "test_add_high_within_cv_cohorts_summary.txt"), sep = "\t")
fwrite(per_study_expanded, file.path(out_dir, "test_add_high_within_cv_cohorts_expanded_per_study_auc.txt"), sep = "\t")
fwrite(pred_expanded, file.path(out_dir, "test_add_high_within_cv_cohorts_expanded_predictions.txt"), sep = "\t")

print(summary_dt)
print(per_study_expanded)
