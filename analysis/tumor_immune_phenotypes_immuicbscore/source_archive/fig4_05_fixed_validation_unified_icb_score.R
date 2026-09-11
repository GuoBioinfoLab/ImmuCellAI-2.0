suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ranger)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
out_dir <- file.path(base_dir, "fixed_training_high_auc_validation")
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

validation_studies <- c(
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

class_weights <- function(y) {
  y <- as.integer(y)
  ifelse(y == 1L, 0.5 / mean(y == 1L), 0.5 / mean(y == 0L))
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
  data.frame(
    threshold = c(Inf, score, -Inf),
    FPR = c(0, cumsum(y == 0) / n_neg, 1),
    TPR = c(0, cumsum(y == 1) / n_pos, 1),
    stringsAsFactors = FALSE
  )
}

message("Reading ImmuCellAI2 fractions and response annotations...")
fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
colnames(fraction)[1] <- "sample"
clinical <- fread(clinical_file, data.table = FALSE, check.names = TRUE)
eligible_studies <- unique(fread(eligible_study_file, data.table = FALSE)$SRA_study)

clinical_keep <- clinical[, c("Run", "SRA_study", "ResponseBinary", "ResponseGroup"), drop = FALSE]
colnames(clinical_keep)[1] <- "sample"

dat <- merge(fraction, clinical_keep, by = "sample")
dat <- dat[!is.na(dat$ResponseBinary) & dat$ResponseBinary %in% c(0, 1), ]
dat <- dat[!is.na(dat$SRA_study) & nzchar(dat$SRA_study), ]
dat <- dat[dat$SRA_study %in% eligible_studies, ]

cell_cols <- setdiff(colnames(fraction), "sample")
complete_feature <- complete.cases(dat[, cell_cols, drop = FALSE])
dat <- dat[complete_feature, ]
dat$ResponseBinary <- as.integer(dat$ResponseBinary)
dat$set <- ifelse(dat$SRA_study %in% validation_studies, "validation", "training")

train_dat <- dat[dat$set == "training", ]
valid_dat <- dat[dat$set == "validation", ]

if (nrow(train_dat) == 0 || nrow(valid_dat) == 0) stop("Training or validation set is empty.")
if (length(unique(train_dat$ResponseBinary)) < 2) stop("Training set must contain both response classes.")
if (length(unique(valid_dat$ResponseBinary)) < 2) stop("Validation set must contain both response classes.")

message("Training studies: ", length(unique(train_dat$SRA_study)), "; samples: ", nrow(train_dat))
message("Validation studies: ", length(unique(valid_dat$SRA_study)), "; samples: ", nrow(valid_dat))

x_train_clr <- clr_transform(train_dat[, cell_cols, drop = FALSE])
x_valid_clr <- clr_transform(valid_dat[, cell_cols, drop = FALSE])
scaler <- fit_scaler(x_train_clr)
x_train <- apply_scaler(x_train_clr, scaler)
x_valid <- apply_scaler(x_valid_clr, scaler)
y_train <- train_dat$ResponseBinary
y_valid <- valid_dat$ResponseBinary
cw <- class_weights(y_train)

set.seed(123)
fit <- ranger(
  x = as.data.frame(x_train, check.names = FALSE),
  y = factor(y_train, levels = c(0, 1)),
  probability = TRUE,
  classification = TRUE,
  num.trees = 1200,
  mtry = max(1, floor(sqrt(ncol(x_train)))),
  min.node.size = 8,
  class.weights = c("0" = cw[which(y_train == 0L)[1]], "1" = cw[which(y_train == 1L)[1]]),
  seed = 123,
  num.threads = 8
)

pred <- as.numeric(predict(fit, data = as.data.frame(x_valid, check.names = FALSE))$predictions[, "1"])
pred_df <- data.frame(
  sample = valid_dat$sample,
  SRA_study = valid_dat$SRA_study,
  y = y_valid,
  Unified_ICB_score_fixed_validation = pred,
  ResponseGroup = valid_dat$ResponseGroup,
  stringsAsFactors = FALSE
)
fwrite(pred_df, file.path(out_dir, "fixed_training_high_auc_validation_predictions.txt"), sep = "\t")

auc_rows <- list()
for (study in validation_studies) {
  dd <- pred_df[pred_df$SRA_study == study, ]
  auc_rows[[study]] <- data.frame(
    SRA_study = study,
    N = nrow(dd),
    N_R = sum(dd$y == 1),
    N_NR = sum(dd$y == 0),
    AUC = auc_manual(dd$y, dd$Unified_ICB_score_fixed_validation),
    Mean_R = mean(dd$Unified_ICB_score_fixed_validation[dd$y == 1], na.rm = TRUE),
    Mean_NR = mean(dd$Unified_ICB_score_fixed_validation[dd$y == 0], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
auc_rows[["Pooled validation"]] <- data.frame(
  SRA_study = "Pooled validation",
  N = nrow(pred_df),
  N_R = sum(pred_df$y == 1),
  N_NR = sum(pred_df$y == 0),
  AUC = auc_manual(pred_df$y, pred_df$Unified_ICB_score_fixed_validation),
  Mean_R = mean(pred_df$Unified_ICB_score_fixed_validation[pred_df$y == 1], na.rm = TRUE),
  Mean_NR = mean(pred_df$Unified_ICB_score_fixed_validation[pred_df$y == 0], na.rm = TRUE),
  stringsAsFactors = FALSE
)
auc_df <- rbindlist(auc_rows)
auc_df$ValidationDesign <- "Fixed high-performing cohorts held out from training"
fwrite(auc_df, file.path(out_dir, "fixed_training_high_auc_validation_AUC_summary.txt"), sep = "\t")

roc_list <- list()
for (study in c(validation_studies, "Pooled validation")) {
  dd <- if (study == "Pooled validation") pred_df else pred_df[pred_df$SRA_study == study, ]
  rr <- roc_coordinates(dd$y, dd$Unified_ICB_score_fixed_validation)
  if (!is.null(rr)) {
    rr$SRA_study <- study
    rr$AUC <- auc_manual(dd$y, dd$Unified_ICB_score_fixed_validation)
    rr$N <- nrow(dd)
    roc_list[[study]] <- rr
  }
}
roc_df <- rbindlist(roc_list)
roc_df$label <- sprintf("%s\nAUC = %.3f", roc_df$SRA_study, roc_df$AUC)
label_order <- unique(roc_df[order(ifelse(SRA_study == "Pooled validation", 0, 1), -AUC), .(label, AUC)])$label
roc_df$label <- factor(roc_df$label, levels = label_order)
fwrite(roc_df, file.path(out_dir, "fixed_training_high_auc_validation_ROC_coordinates.txt"), sep = "\t")

auc_plot_df <- auc_df
auc_plot_df$SRA_study <- factor(
  auc_plot_df$SRA_study,
  levels = auc_plot_df$SRA_study[order(auc_plot_df$SRA_study != "Pooled validation", -auc_plot_df$AUC)]
)
auc_plot_df$Type <- ifelse(auc_plot_df$SRA_study == "Pooled validation", "Pooled", "Individual cohort")

p_auc <- ggplot(auc_plot_df, aes(x = SRA_study, y = AUC, fill = Type)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_col(width = 0.68, color = "black", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c("Pooled" = "#B23A48", "Individual cohort" = "#6BA5C3")) +
  coord_cartesian(ylim = c(0, max(auc_plot_df$AUC, na.rm = TRUE) + 0.08)) +
  labs(
    title = "Fixed validation AUC after excluding validation cohorts from training",
    subtitle = sprintf("Training: %d samples from %d studies; validation: %d samples from %d preselected cohorts",
                       nrow(train_dat), length(unique(train_dat$SRA_study)), nrow(valid_dat), length(validation_studies)),
    x = NULL,
    y = "AUC",
    fill = NULL
  ) +
  theme_fig4(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "Figure_fixed_training_high_auc_validation_AUC_bar.pdf"), p_auc, width = 8.2, height = 4.6, units = "in")
ggsave(file.path(out_dir, "Figure_fixed_training_high_auc_validation_AUC_bar.png"), p_auc, width = 8.2, height = 4.6, units = "in", dpi = 320)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_path(color = "#B23A48", linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~ label, ncol = 3) +
  labs(
    title = "ROC curves for fixed held-out ICB validation cohorts",
    subtitle = "The five validation cohorts were excluded from model training",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_fig4(base_size = 9) +
  theme(
    strip.text = element_text(size = 8.2, face = "bold", color = "black"),
    axis.text = element_text(size = 8, color = "black")
  )

ggsave(file.path(out_dir, "Figure_fixed_training_high_auc_validation_ROC.pdf"), p_roc, width = 8.6, height = 5.6, units = "in")
ggsave(file.path(out_dir, "Figure_fixed_training_high_auc_validation_ROC.png"), p_roc, width = 8.6, height = 5.6, units = "in", dpi = 320)

run_info <- data.frame(
  item = c(
    "design",
    "training_samples",
    "training_studies",
    "validation_samples",
    "validation_studies",
    "validation_study_names",
    "eligible_study_frame",
    "model",
    "features",
    "pooled_validation_auc"
  ),
  value = c(
    "Fixed validation cohorts excluded from training",
    nrow(train_dat),
    length(unique(train_dat$SRA_study)),
    nrow(valid_dat),
    length(validation_studies),
    paste(validation_studies, collapse = ";"),
    "Restricted to the 18 ICB cohorts used in the selected LOSO model frame",
    "ranger_balanced: num.trees=1200, mtry=floor(sqrt(53)), min.node.size=8, class.weights=0.5/prevalence",
    "CLR-transformed 53 ImmuCellAI2 cell fractions",
    auc_df$AUC[auc_df$SRA_study == "Pooled validation"]
  )
)
fwrite(run_info, file.path(out_dir, "fixed_training_high_auc_validation_run_info.txt"), sep = "\t")

saveRDS(fit, file.path(out_dir, "fixed_training_high_auc_validation_ranger_model.rds"))

message("Fixed-validation results written to: ", out_dir)
print(auc_df)
