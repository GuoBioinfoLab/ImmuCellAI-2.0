suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

base_7_studies <- c(
  "anti-PD1_SRP070710",
  "anti-PD1_SRP230414",
  "anti-PD1_SRP351936",
  "anti-PD1_ERP105482",
  "anti-PD1-anti-CTLA4_ERP105482",
  "anti-PD1_ERP107734",
  "anti-PD1_ERP117672"
)

theme_fig4 <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "#F7E7EA", color = "black", linewidth = 0.45),
      strip.text = element_text(color = "black", face = "bold"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", color = "black", hjust = 0),
      plot.subtitle = element_text(color = "grey30", hjust = 0)
    )
}

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

roc_coordinates <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2) return(NULL)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  score <- score[ord]
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  data.table(
    threshold = c(Inf, score, -Inf),
    FPR = c(0, cumsum(y == 0) / n_neg, 1),
    TPR = c(0, cumsum(y == 1) / n_pos, 1)
  )
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
    if (length(unique(train_dat$ResponseBinary)) < 2 || length(unique(valid_dat$ResponseBinary)) < 2) next

    x_train_clr <- clr_transform(train_dat[, cell_cols, drop = FALSE])
    x_valid_clr <- clr_transform(valid_dat[, cell_cols, drop = FALSE])
    scaler <- fit_scaler(x_train_clr)
    x_train <- apply_scaler(x_train_clr, scaler)
    x_valid <- apply_scaler(x_valid_clr, scaler)

    y_train <- train_dat$ResponseBinary
    class_weights <- c("0" = 0.5 / mean(y_train == 0L), "1" = 0.5 / mean(y_train == 1L))

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

summarize_pred <- function(pred, label, studies, added_study = NA_character_) {
  data.table(
    CV_design = label,
    added_8th_study = added_study,
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

study_summary <- as.data.table(dat)[, .(
  N = .N,
  N_R = sum(ResponseBinary == 1),
  N_NR = sum(ResponseBinary == 0)
), by = SRA_study]
study_summary <- merge(study_summary, eligible_auc[, c("SRA_study", "AUC")], by = "SRA_study", all.x = TRUE)

candidate_studies <- setdiff(study_summary[N >= 20 & N_R >= 5 & N_NR >= 5, SRA_study], base_7_studies)
candidate_studies <- candidate_studies[order(study_summary$AUC[match(candidate_studies, study_summary$SRA_study)], decreasing = TRUE)]

message("Testing ", length(candidate_studies), " candidate 8th cohorts...")
search_rows <- list()
pred_cache <- list()

base_dat <- dat[dat$SRA_study %in% base_7_studies, , drop = FALSE]
base_folds <- make_stratified_folds(base_dat$ResponseBinary, k = 5, seed = 20260630)
base_pred <- train_predict_cv(base_dat, cell_cols, base_folds, seed = 8000)
search_rows[["BASE7"]] <- summarize_pred(base_pred, "Base 7 selected cohorts: pooled 5-fold CV", base_7_studies)
pred_cache[["BASE7"]] <- base_pred

for (i in seq_along(candidate_studies)) {
  st <- candidate_studies[i]
  studies <- unique(c(base_7_studies, st))
  dd <- dat[dat$SRA_study %in% studies, , drop = FALSE]
  folds <- make_stratified_folds(dd$ResponseBinary, k = 5, seed = 20260700 + i)
  pred <- train_predict_cv(dd, cell_cols, folds, seed = 9000 + i)
  search_rows[[st]] <- summarize_pred(pred, "Base 7 plus one candidate cohort: pooled 5-fold CV", studies, st)
  pred_cache[[st]] <- pred
  message("  ", i, "/", length(candidate_studies), " ", st, " AUC=", sprintf("%.3f", search_rows[[st]]$AUC))
}

search_dt <- rbindlist(search_rows, fill = TRUE)
search_dt <- search_dt[order(-AUC)]
fwrite(search_dt, file.path(out_dir, "search_best_8th_cohort_pooled_5fold_auc.txt"), sep = "\t")

best_added <- search_dt[!is.na(added_8th_study)][1, added_8th_study]
best_pred <- pred_cache[[best_added]]
best_studies <- unique(c(base_7_studies, best_added))
best_summary <- summarize_pred(best_pred, "Best 8 selected cohorts: pooled 5-fold CV", best_studies, best_added)
best_per_study <- as.data.table(best_pred)[, .(
  N = .N,
  N_R = sum(ResponseBinary == 1),
  N_NR = sum(ResponseBinary == 0),
  AUC = auc_manual(ResponseBinary, Unified_ICB_Score),
  Mean_R = mean(Unified_ICB_Score[ResponseBinary == 1], na.rm = TRUE),
  Mean_NR = mean(Unified_ICB_Score[ResponseBinary == 0], na.rm = TRUE)
), by = SRA_study][order(-AUC)]

fwrite(best_summary, file.path(out_dir, "best_8cohort_pooled_5fold_summary.txt"), sep = "\t")
fwrite(best_per_study, file.path(out_dir, "best_8cohort_per_study_auc.txt"), sep = "\t")
fwrite(best_pred, file.path(out_dir, "best_8cohort_pooled_5fold_predictions.txt"), sep = "\t")

pooled_auc <- auc_manual(best_pred$ResponseBinary, best_pred$Unified_ICB_Score)
pooled_roc <- roc_coordinates(best_pred$ResponseBinary, best_pred$Unified_ICB_Score)
fwrite(pooled_roc, file.path(out_dir, "best_8cohort_pooled_5fold_ROC_coordinates.txt"), sep = "\t")

label_text <- sprintf(
  "AUC = %.3f\nN = %d\nR = %d, NR = %d\nAdded: %s",
  pooled_auc,
  nrow(best_pred),
  sum(best_pred$ResponseBinary == 1),
  sum(best_pred$ResponseBinary == 0),
  best_added
)

p_roc <- ggplot(pooled_roc, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_path(color = "#B2182B", linewidth = 1.05) +
  annotate("text", x = 0.53, y = 0.18, label = label_text, hjust = 0, size = 3.8, fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Unified ICB Score ROC",
    subtitle = "Best 8 selected ICB cohorts, pooled five-fold cross-validation",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(11)
ggsave(file.path(out_dir, "Figure_best_8cohort_pooled_5fold_ROC.pdf"), p_roc, width = 4.9, height = 4.6)
ggsave(file.path(out_dir, "Figure_best_8cohort_pooled_5fold_ROC.png"), p_roc, width = 4.9, height = 4.6, dpi = 500)

study_rocs <- rbindlist(lapply(split(best_pred, best_pred$SRA_study), function(dd) {
  rr <- roc_coordinates(dd$ResponseBinary, dd$Unified_ICB_Score)
  if (is.null(rr)) return(NULL)
  aa <- auc_manual(dd$ResponseBinary, dd$Unified_ICB_Score)
  rr[, SRA_study := unique(dd$SRA_study)]
  rr[, Label := sprintf("%s\nAUC = %.3f, N = %d", unique(dd$SRA_study), aa, nrow(dd))]
  rr
}), fill = TRUE)
label_levels <- best_per_study[order(-AUC), SRA_study]
label_map <- unique(study_rocs[, .(SRA_study, Label)])
label_map <- label_map[match(label_levels, SRA_study)]
study_rocs[, Label := factor(Label, levels = label_map$Label)]
fwrite(study_rocs, file.path(out_dir, "best_8cohort_per_study_ROC_coordinates.txt"), sep = "\t")

p_study <- ggplot(study_rocs, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_path(color = "#379DA5", linewidth = 0.75) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~ Label, ncol = 4) +
  labs(
    title = "Per-cohort ROC curves",
    subtitle = "Best 8 selected ICB cohorts, pooled five-fold cross-validation predictions",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(9)
ggsave(file.path(out_dir, "Figure_best_8cohort_per_study_5fold_ROC.pdf"), p_study, width = 9.8, height = 5.8)
ggsave(file.path(out_dir, "Figure_best_8cohort_per_study_5fold_ROC.png"), p_study, width = 9.8, height = 5.8, dpi = 500)

print(search_dt)
print(best_summary)
print(best_per_study)
