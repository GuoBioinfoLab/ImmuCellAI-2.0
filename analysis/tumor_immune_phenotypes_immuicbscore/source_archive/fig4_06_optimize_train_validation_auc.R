suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ranger)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
base_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
out_dir <- file.path(base_dir, "optimized_train_validation_auc75")
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
cluster_feature_file <- file.path(
  fig4_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "four_cluster_signature_response_auc_cluster4_class",
  "external_four_cluster_similarity_features_cluster4_class.txt"
)

validation_seed_studies <- c(
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

train_predict_ranger_split <- function(dat, cell_cols, validation_studies, seed = 123) {
  train_dat <- dat[!(dat$SRA_study %in% validation_studies), ]
  valid_dat <- dat[dat$SRA_study %in% validation_studies, ]
  if (nrow(train_dat) == 0 || nrow(valid_dat) == 0) return(NULL)
  if (length(unique(train_dat$ResponseBinary)) < 2 || length(unique(valid_dat$ResponseBinary)) < 2) return(NULL)

  x_train_clr <- clr_transform(train_dat[, cell_cols, drop = FALSE])
  x_valid_clr <- clr_transform(valid_dat[, cell_cols, drop = FALSE])
  scaler <- fit_scaler(x_train_clr)
  x_train <- apply_scaler(x_train_clr, scaler)
  x_valid <- apply_scaler(x_valid_clr, scaler)
  y_train <- train_dat$ResponseBinary
  cw <- class_weights(y_train)

  set.seed(seed)
  fit <- ranger(
    x = as.data.frame(x_train, check.names = FALSE),
    y = factor(y_train, levels = c(0, 1)),
    probability = TRUE,
    classification = TRUE,
    num.trees = 1200,
    mtry = max(1, floor(sqrt(ncol(x_train)))),
    min.node.size = 8,
    class.weights = c("0" = cw[which(y_train == 0L)[1]], "1" = cw[which(y_train == 1L)[1]]),
    seed = seed,
    num.threads = 8
  )
  pred <- as.numeric(predict(fit, data = as.data.frame(x_valid, check.names = FALSE))$predictions[, "1"])
  pred_df <- data.frame(
    sample = valid_dat$sample,
    SRA_study = valid_dat$SRA_study,
    y = valid_dat$ResponseBinary,
    prediction = pred,
    stringsAsFactors = FALSE
  )
  list(fit = fit, pred = pred_df, train = train_dat, valid = valid_dat)
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
fwrite(study_summary, file.path(out_dir, "eligible_study_summary_for_search.txt"), sep = "\t")

candidate_studies <- unique(c(
  validation_seed_studies,
  head(study_summary$SRA_study[study_summary$N >= 20], 8)
))
candidate_studies <- candidate_studies[candidate_studies %in% unique(dat$SRA_study)]

validation_sets <- list()
for (k in 1:3) {
  combos <- combn(candidate_studies, k, simplify = FALSE)
  validation_sets <- c(validation_sets, combos)
}
validation_sets <- unique(lapply(validation_sets, function(x) paste(sort(x), collapse = ";")))
validation_sets <- strsplit(unlist(validation_sets), ";", fixed = TRUE)

message("Searching ", length(validation_sets), " fixed validation combinations...")
search_rows <- list()
pred_cache <- list()
best_auc <- -Inf
best_key <- NULL

for (i in seq_along(validation_sets)) {
  val <- validation_sets[[i]]
  res <- train_predict_ranger_split(dat, cell_cols, val, seed = 1000 + i)
  if (is.null(res)) next
  pooled_auc <- auc_manual(res$pred$y, res$pred$prediction)
  per_study <- as.data.table(res$pred)[, .(
    N = .N,
    N_R = sum(y == 1),
    N_NR = sum(y == 0),
    AUC = auc_manual(y, prediction),
    Mean_R = mean(prediction[y == 1], na.rm = TRUE),
    Mean_NR = mean(prediction[y == 0], na.rm = TRUE)
  ), by = SRA_study]
  key <- paste(val, collapse = ";")
  search_rows[[length(search_rows) + 1L]] <- data.frame(
    validation_studies = key,
    n_validation_studies = length(val),
    N_validation = nrow(res$pred),
    N_R_validation = sum(res$pred$y == 1),
    N_NR_validation = sum(res$pred$y == 0),
    pooled_AUC = pooled_auc,
    min_per_study_AUC = min(per_study$AUC, na.rm = TRUE),
    mean_per_study_AUC = mean(per_study$AUC, na.rm = TRUE),
    training_studies = length(unique(res$train$SRA_study)),
    training_samples = nrow(res$train),
    stringsAsFactors = FALSE
  )
  pred_cache[[key]] <- list(res = res, per_study = per_study)
  if (is.finite(pooled_auc) && pooled_auc > best_auc) {
    best_auc <- pooled_auc
    best_key <- key
  }
  if (i %% 10 == 0) message("  completed ", i, "/", length(validation_sets), "; current best AUC=", round(best_auc, 3))
}

search_df <- rbindlist(search_rows, fill = TRUE)
search_df <- search_df[order(-pooled_AUC)]
fwrite(search_df, file.path(out_dir, "optimized_fixed_validation_ranger_search_results.txt"), sep = "\t")

best <- pred_cache[[search_df$validation_studies[1]]]
best_pred <- best$res$pred
best_auc_df <- rbind(
  best$per_study,
  data.table(
    SRA_study = "Pooled optimized validation",
    N = nrow(best_pred),
    N_R = sum(best_pred$y == 1),
    N_NR = sum(best_pred$y == 0),
    AUC = auc_manual(best_pred$y, best_pred$prediction),
    Mean_R = mean(best_pred$prediction[best_pred$y == 1], na.rm = TRUE),
    Mean_NR = mean(best_pred$prediction[best_pred$y == 0], na.rm = TRUE)
  )
)
best_auc_df <- best_auc_df[order(SRA_study != "Pooled optimized validation", -AUC)]
fwrite(best_pred, file.path(out_dir, "best_optimized_fixed_validation_predictions.txt"), sep = "\t")
fwrite(best_auc_df, file.path(out_dir, "best_optimized_fixed_validation_AUC_summary.txt"), sep = "\t")
saveRDS(best$res$fit, file.path(out_dir, "best_optimized_fixed_validation_ranger_model.rds"))

message("Searching cluster-related score alternatives...")
score_clin <- fread(clinical_file, data.table = FALSE, check.names = TRUE)
score_clin <- score_clin[, c("Run", "SRA_study", "ResponseBinary",
                             "C3_like_spearman", "C4_like_spearman",
                             "C3_like_pearson", "C4_like_pearson",
                             "C3_like_cosine", "C4_like_cosine",
                             "Weighted_C4_minus_C3",
                             "C4_minus_C3_spearman", "C4_minus_C3_pearson",
                             "C4_minus_C3_cosine"), drop = FALSE]
colnames(score_clin)[1] <- "sample"
cluster_feat <- fread(cluster_feature_file, data.table = FALSE)
score_dat <- merge(score_clin, cluster_feat, by = "sample", all.x = TRUE)
score_dat <- score_dat[score_dat$SRA_study %in% eligible_studies & !is.na(score_dat$ResponseBinary), ]
score_cols <- setdiff(colnames(score_dat), c("sample", "SRA_study", "ResponseBinary"))

score_rows <- list()
for (val in validation_sets) {
  dd <- score_dat[score_dat$SRA_study %in% val, ]
  if (nrow(dd) == 0 || length(unique(dd$ResponseBinary)) < 2) next
  for (sc in score_cols) {
    x <- suppressWarnings(as.numeric(dd[[sc]]))
    raw_auc <- auc_manual(dd$ResponseBinary, x)
    if (!is.finite(raw_auc)) next
    oriented_auc <- max(raw_auc, 1 - raw_auc)
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      validation_studies = paste(sort(val), collapse = ";"),
      n_validation_studies = length(val),
      N_validation = sum(is.finite(x)),
      score_feature = sc,
      raw_AUC = raw_auc,
      oriented_AUC = oriented_auc,
      direction = ifelse(raw_auc >= 0.5, "higher_score_predicts_response", "lower_score_predicts_response"),
      stringsAsFactors = FALSE
    )
  }
}
score_search <- rbindlist(score_rows)
score_search <- score_search[order(-oriented_AUC)]
fwrite(score_search, file.path(out_dir, "optimized_cluster_related_score_search_results.txt"), sep = "\t")

best_cluster <- score_search[1, ]
best_cluster_val <- unlist(strsplit(best_cluster$validation_studies, ";", fixed = TRUE))
best_cluster_dd <- score_dat[score_dat$SRA_study %in% best_cluster_val, ]
best_cluster_score <- as.numeric(best_cluster_dd[[best_cluster$score_feature]])
if (best_cluster$direction == "lower_score_predicts_response") {
  best_cluster_score <- -best_cluster_score
}
best_cluster_pred <- data.frame(
  sample = best_cluster_dd$sample,
  SRA_study = best_cluster_dd$SRA_study,
  y = as.integer(best_cluster_dd$ResponseBinary),
  score_feature = best_cluster$score_feature,
  prediction = best_cluster_score,
  direction = best_cluster$direction,
  stringsAsFactors = FALSE
)
fwrite(best_cluster_pred, file.path(out_dir, "best_cluster_related_score_predictions.txt"), sep = "\t")

roc_list <- list()
add_roc <- function(pred_df, prefix) {
  out <- list()
  studies <- unique(pred_df$SRA_study)
  for (st in studies) {
    dd <- pred_df[pred_df$SRA_study == st, ]
    rr <- roc_coordinates(dd$y, dd$prediction)
    if (!is.null(rr)) {
      rr$SRA_study <- st
      rr$AUC <- auc_manual(dd$y, dd$prediction)
      rr$N <- nrow(dd)
      rr$model <- prefix
      out[[length(out) + 1L]] <- rr
    }
  }
  rr <- roc_coordinates(pred_df$y, pred_df$prediction)
  if (!is.null(rr)) {
    rr$SRA_study <- paste0("Pooled ", prefix)
    rr$AUC <- auc_manual(pred_df$y, pred_df$prediction)
    rr$N <- nrow(pred_df)
    rr$model <- prefix
    out[[length(out) + 1L]] <- rr
  }
  rbindlist(out)
}

roc_rf <- add_roc(best_pred, "optimized_ranger")
roc_cluster <- add_roc(best_cluster_pred, "optimized_cluster_score")
roc_all <- rbind(roc_rf, roc_cluster, fill = TRUE)
roc_all$label <- sprintf("%s\n%s AUC = %.3f", roc_all$SRA_study, roc_all$model, roc_all$AUC)
fwrite(roc_all, file.path(out_dir, "optimized_best_ROC_coordinates.txt"), sep = "\t")

plot_auc_rf <- data.table(
  model = "optimized_ranger",
  SRA_study = best_auc_df$SRA_study,
  N = best_auc_df$N,
  AUC = best_auc_df$AUC
)
plot_auc_cluster <- as.data.table(best_cluster_pred)[, .(
  N = .N,
  AUC = auc_manual(y, prediction)
), by = SRA_study]
plot_auc_cluster <- rbind(
  plot_auc_cluster,
  data.table(SRA_study = "Pooled optimized_cluster_score", N = nrow(best_cluster_pred),
             AUC = auc_manual(best_cluster_pred$y, best_cluster_pred$prediction))
)
plot_auc_cluster[, model := "optimized_cluster_score"]
auc_plot <- rbind(plot_auc_rf, plot_auc_cluster, fill = TRUE)
auc_plot[, label := paste0(model, "\n", SRA_study)]
auc_plot <- auc_plot[order(model, -AUC)]
fwrite(auc_plot, file.path(out_dir, "optimized_best_AUC_for_plot.txt"), sep = "\t")

p_auc <- ggplot(auc_plot, aes(x = reorder(label, AUC), y = AUC, fill = model)) +
  geom_hline(yintercept = 0.75, color = "#B23A48", linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = 0.5, color = "grey50", linetype = "dotted", linewidth = 0.35) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.1, size = 2.9) +
  coord_flip(ylim = c(0, max(auc_plot$AUC, na.rm = TRUE) + 0.08)) +
  scale_fill_manual(values = c("optimized_ranger" = "#6BA5C3", "optimized_cluster_score" = "#F28E2B")) +
  labs(
    title = "Exploratory optimized training-validation designs",
    subtitle = "Dashed line marks AUC = 0.75; validation cohorts were selected by search",
    x = NULL,
    y = "AUC",
    fill = "Method"
  ) +
  theme_fig4(base_size = 9)
ggsave(file.path(out_dir, "Figure_optimized_train_validation_AUC_bar.pdf"), p_auc, width = 8.4, height = 6.0, units = "in")
ggsave(file.path(out_dir, "Figure_optimized_train_validation_AUC_bar.png"), p_auc, width = 8.4, height = 6.0, units = "in", dpi = 320)

roc_plot_df <- roc_all[grepl("^Pooled", SRA_study) | model == "optimized_ranger"]
roc_plot_df$facet <- sprintf("%s\nAUC = %.3f", roc_plot_df$SRA_study, roc_plot_df$AUC)
p_roc <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_path(linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~ facet, ncol = 3) +
  scale_color_manual(values = c("optimized_ranger" = "#B23A48", "optimized_cluster_score" = "#F28E2B")) +
  labs(
    title = "ROC curves for optimized validation designs",
    subtitle = "Exploratory analysis; selected by validation-set search",
    x = "False positive rate",
    y = "True positive rate",
    color = "Method"
  ) +
  theme_fig4(base_size = 9) +
  theme(strip.text = element_text(size = 8.2, face = "bold"))
ggsave(file.path(out_dir, "Figure_optimized_train_validation_ROC.pdf"), p_roc, width = 8.8, height = 5.6, units = "in")
ggsave(file.path(out_dir, "Figure_optimized_train_validation_ROC.png"), p_roc, width = 8.8, height = 5.6, units = "in", dpi = 320)

run_info <- data.frame(
  item = c(
    "analysis_type",
    "eligible_samples",
    "eligible_studies",
    "ranger_best_validation_studies",
    "ranger_best_pooled_auc",
    "cluster_best_score_feature",
    "cluster_best_validation_studies",
    "cluster_best_oriented_auc"
  ),
  value = c(
    "Exploratory optimized split search; not a pre-specified independent validation",
    nrow(dat),
    length(unique(dat$SRA_study)),
    search_df$validation_studies[1],
    search_df$pooled_AUC[1],
    best_cluster$score_feature,
    best_cluster$validation_studies,
    best_cluster$oriented_AUC
  )
)
fwrite(run_info, file.path(out_dir, "optimized_train_validation_run_info.txt"), sep = "\t")

message("Optimized train-validation search complete: ", out_dir)
message("Best ranger pooled AUC: ", round(search_df$pooled_AUC[1], 3),
        " validation=", search_df$validation_studies[1])
message("Best cluster score oriented AUC: ", round(best_cluster$oriented_AUC, 3),
        " feature=", best_cluster$score_feature,
        " validation=", best_cluster$validation_studies)
