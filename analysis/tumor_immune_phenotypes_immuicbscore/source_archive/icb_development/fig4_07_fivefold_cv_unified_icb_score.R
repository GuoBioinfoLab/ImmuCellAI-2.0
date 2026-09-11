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

high_loso_studies <- c(
  "anti-PD1_SRP070710",
  "anti-PD1_SRP230414",
  "anti-PD1_SRP351936",
  "anti-PD1_ERP105482",
  "anti-PD1-anti-CTLA4_ERP105482"
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
    if (length(unique(train_dat$ResponseBinary)) < 2 || length(unique(valid_dat$ResponseBinary)) < 2) {
      next
    }

    x_train_clr <- clr_transform(train_dat[, cell_cols, drop = FALSE])
    x_valid_clr <- clr_transform(valid_dat[, cell_cols, drop = FALSE])
    scaler <- fit_scaler(x_train_clr)
    x_train <- apply_scaler(x_train_clr, scaler)
    x_valid <- apply_scaler(x_valid_clr, scaler)
    y_train <- train_dat$ResponseBinary
    prevalence_pos <- mean(y_train == 1L)
    prevalence_neg <- mean(y_train == 0L)
    class_weights <- c("0" = 0.5 / prevalence_neg, "1" = 0.5 / prevalence_pos)

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

summarize_auc <- function(pred, label) {
  data.table(
    CV_design = label,
    N = nrow(pred),
    N_R = sum(pred$ResponseBinary == 1),
    N_NR = sum(pred$ResponseBinary == 0),
    AUC = auc_manual(pred$ResponseBinary, pred$Unified_ICB_Score),
    Mean_R = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 1], na.rm = TRUE),
    Mean_NR = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 0], na.rm = TRUE)
  )
}

message("Reading external ICB data...")
fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
colnames(fraction)[1] <- "sample"
clinical <- fread(clinical_file, data.table = FALSE, check.names = TRUE)
eligible_auc <- fread(eligible_study_file, data.table = FALSE)
eligible_studies <- unique(eligible_auc$SRA_study)

clinical_keep <- clinical[, c("Run", "SRA_study", "ResponseBinary", "ResponseGroup"), drop = FALSE]
colnames(clinical_keep)[1] <- "sample"
dat <- merge(fraction, clinical_keep, by = "sample")
dat <- dat[!is.na(dat$ResponseBinary) & dat$ResponseBinary %in% c(0, 1), ]
dat <- dat[!is.na(dat$SRA_study) & dat$SRA_study %in% eligible_studies, ]
cell_cols <- setdiff(colnames(fraction), "sample")
dat <- dat[complete.cases(dat[, cell_cols, drop = FALSE]), ]
dat$ResponseBinary <- as.integer(dat$ResponseBinary)

study_summary <- as.data.table(dat)[, .(
  N = .N,
  N_R = sum(ResponseBinary == 1),
  N_NR = sum(ResponseBinary == 0)
), by = SRA_study]
study_summary <- merge(study_summary, eligible_auc[, c("SRA_study", "AUC")], by = "SRA_study", all.x = TRUE)
study_summary <- study_summary[order(-AUC)]
fwrite(study_summary, file.path(out_dir, "eligible_study_summary_for_fivefold_cv.txt"), sep = "\t")

message("Running pooled 5-fold CV on all eligible ICB cohorts...")
folds_all <- make_stratified_folds(dat$ResponseBinary, k = 5, seed = 20260628)
pred_all <- train_predict_cv(dat, cell_cols, folds_all, seed = 3000, num_trees = 1200)
fwrite(pred_all, file.path(out_dir, "pooled_all_18studies_5fold_predictions.txt"), sep = "\t")

message("Running pooled 5-fold CV on high-LOSO cohorts...")
dat_high <- dat[dat$SRA_study %in% high_loso_studies, , drop = FALSE]
folds_high <- make_stratified_folds(dat_high$ResponseBinary, k = 5, seed = 20260629)
pred_high <- train_predict_cv(dat_high, cell_cols, folds_high, seed = 4000, num_trees = 1200)
fwrite(pred_high, file.path(out_dir, "pooled_high_LOSO_5studies_5fold_predictions.txt"), sep = "\t")

message("Running within-study 5-fold CV for eligible individual cohorts...")
within_rows <- list()
within_pred_rows <- list()
study_order <- study_summary$SRA_study
for (i in seq_along(study_order)) {
  st <- study_order[i]
  dd <- dat[dat$SRA_study == st, , drop = FALSE]
  if (nrow(dd) < 25 || min(table(dd$ResponseBinary)) < 5) next
  folds <- make_stratified_folds(dd$ResponseBinary, k = 5, seed = 5000 + i)
  pred <- train_predict_cv(dd, cell_cols, folds, seed = 6000 + i, num_trees = 1200)
  if (nrow(pred) == 0) next
  within_pred_rows[[st]] <- pred
  within_rows[[st]] <- data.table(
    SRA_study = st,
    N = nrow(pred),
    N_R = sum(pred$ResponseBinary == 1),
    N_NR = sum(pred$ResponseBinary == 0),
    AUC = auc_manual(pred$ResponseBinary, pred$Unified_ICB_Score),
    Mean_R = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 1], na.rm = TRUE),
    Mean_NR = mean(pred$Unified_ICB_Score[pred$ResponseBinary == 0], na.rm = TRUE),
    LOSO_AUC = study_summary$AUC[match(st, study_summary$SRA_study)]
  )
}
within_auc <- rbindlist(within_rows, fill = TRUE)
within_auc <- within_auc[order(-AUC)]
within_pred <- rbindlist(within_pred_rows, fill = TRUE)
fwrite(within_auc, file.path(out_dir, "within_study_5fold_auc_by_study.txt"), sep = "\t")
fwrite(within_pred, file.path(out_dir, "within_study_5fold_predictions.txt"), sep = "\t")

summary_auc <- rbindlist(list(
  summarize_auc(pred_all, "All eligible 18 ICB studies: pooled 5-fold CV"),
  summarize_auc(pred_high, "High-LOSO 5 ICB studies: pooled 5-fold CV")
), fill = TRUE)
fwrite(summary_auc, file.path(out_dir, "pooled_5fold_auc_summary.txt"), sep = "\t")

above075 <- rbindlist(list(
  data.table(ValidationType = "LOSO external", eligible_auc[eligible_auc$AUC >= 0.75, c("SRA_study", "N", "N_R", "N_NR", "AUC")]),
  data.table(ValidationType = "Within-study 5-fold CV", within_auc[within_auc$AUC >= 0.75, c("SRA_study", "N", "N_R", "N_NR", "AUC")])
), fill = TRUE)
above075 <- above075[order(ValidationType, -AUC)]
fwrite(above075, file.path(out_dir, "studies_with_auc_ge_0.75_LOSO_or_5fold.txt"), sep = "\t")

plot_auc <- rbindlist(list(
  data.table(Design = summary_auc$CV_design, AUC = summary_auc$AUC, N = summary_auc$N),
  data.table(Design = paste0("Within-study CV: ", within_auc$SRA_study), AUC = within_auc$AUC, N = within_auc$N)
), fill = TRUE)
plot_auc$Design <- factor(plot_auc$Design, levels = rev(plot_auc$Design[order(plot_auc$AUC)]))

p_bar <- ggplot(plot_auc, aes(x = Design, y = AUC, fill = AUC >= 0.75)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#B2182B", linewidth = 0.45) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c("FALSE" = "#95A8AC", "TRUE" = "#F66463"), guide = "none") +
  labs(
    title = "Unified ICB Score five-fold cross-validation",
    subtitle = "CLR-transformed 53 ImmuCellAI2 fractions + balanced ranger",
    x = NULL,
    y = "AUC"
  ) +
  theme_fig4(9)
ggsave(file.path(out_dir, "Figure_fivefold_cv_unified_icb_score_auc_bar.pdf"), p_bar, width = 7.2, height = 5.8)
ggsave(file.path(out_dir, "Figure_fivefold_cv_unified_icb_score_auc_bar.png"), p_bar, width = 7.2, height = 5.8, dpi = 450)

writeLines(c(
  "Unified ICB Score 5-fold CV analysis",
  "Method: CLR-transformed 53 ImmuCellAI2 fractions + training-fold scaling + balanced ranger random forest.",
  paste0("Eligible studies: ", length(unique(dat$SRA_study)), "; samples: ", nrow(dat),
         "; responders: ", sum(dat$ResponseBinary == 1), "; non-responders: ", sum(dat$ResponseBinary == 0), "."),
  "Output files:",
  "  pooled_5fold_auc_summary.txt",
  "  within_study_5fold_auc_by_study.txt",
  "  studies_with_auc_ge_0.75_LOSO_or_5fold.txt",
  "  Figure_fivefold_cv_unified_icb_score_auc_bar.pdf/png"
), con = file.path(out_dir, "README_fivefold_cv_unified_icb_score.txt"))

message("Done. Results written to: ", out_dir)
