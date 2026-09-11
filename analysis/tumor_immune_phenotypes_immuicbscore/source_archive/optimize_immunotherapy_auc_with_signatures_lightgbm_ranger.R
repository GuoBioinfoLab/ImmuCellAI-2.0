options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
base_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")

score_file <- file.path(base_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt")
fraction_file <- file.path(base_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
expr_file <- file.path(cluster_dir, "RNAseqTpmSymbol.txt")
out_dir <- file.path(base_dir, "optimized_models_auc_signature_lightgbm_ranger")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(ranger)
  library(lightgbm)
  library(pROC)
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

read_matrix_gene_first <- function(file, select_samples = NULL) {
  con <- file(file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  second <- readLines(con, n = 1L, warn = FALSE)
  first_fields <- strsplit(first, "\t", fixed = TRUE)[[1L]]
  second_fields <- strsplit(second, "\t", fixed = TRUE)[[1L]]

  if (length(first_fields) == length(second_fields) - 1L) {
    keep <- if (is.null(select_samples)) seq_along(first_fields) else match(select_samples, first_fields)
    keep <- keep[!is.na(keep)]
    dt <- fread(
      file,
      sep = "\t",
      header = FALSE,
      skip = 1L,
      select = c(1L, keep + 1L),
      col.names = c("Gene", first_fields[keep]),
      check.names = FALSE,
      data.table = FALSE
    )
  } else {
    dt <- fread(file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
    if (!is.null(select_samples)) {
      cols <- c(colnames(dt)[1], intersect(select_samples, colnames(dt)))
      dt <- dt[, cols, drop = FALSE]
    }
  }

  gene <- as.character(dt[[1L]])
  dt[[1L]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "double"
  rownames(mat) <- toupper(gene)
  mat[!is.finite(mat)] <- 0
  if (anyDuplicated(rownames(mat))) mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  mat
}

clean_response <- function(score_df) {
  y <- rep(NA_integer_, nrow(score_df))
  if ("ResponseGroup" %in% colnames(score_df)) {
    y[score_df$ResponseGroup == "Responder"] <- 1L
    y[score_df$ResponseGroup == "Non-responder"] <- 0L
  }
  if ("Response" %in% colnames(score_df)) {
    y[is.na(y) & score_df$Response == "R"] <- 1L
    y[is.na(y) & score_df$Response == "NR"] <- 0L
  }
  if ("ResponseBinary" %in% colnames(score_df)) {
    rb <- as.character(score_df$ResponseBinary)
    y[is.na(y) & rb == "1"] <- 1L
    y[is.na(y) & rb == "0"] <- 0L
  }
  y
}

one_hot_meta <- function(meta, cols) {
  dat <- meta[, cols, drop = FALSE]
  for (cc in names(dat)) {
    dat[[cc]] <- as.character(dat[[cc]])
    dat[[cc]][is.na(dat[[cc]]) | dat[[cc]] == ""] <- "Unknown"
    top <- names(sort(table(dat[[cc]]), decreasing = TRUE))[seq_len(min(12L, length(unique(dat[[cc]]))))]
    dat[[cc]][!dat[[cc]] %in% top] <- "Other"
    dat[[cc]] <- factor(dat[[cc]])
  }
  mm <- model.matrix(~ . - 1, data = dat)
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  mm
}

signature_sets <- list(
  IFNG_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP5"),
  Cytotoxicity = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "CTSW", "FGFBP2"),
  Tcell_inflamed_GEP = c("CD8A", "CD8B", "CXCL9", "CXCL10", "IDO1", "IFNG", "HLA-DRA", "STAT1", "GZMA", "PRF1", "LAG3", "PDCD1", "TIGIT", "CMKLR1", "CCL5", "CXCR6", "NKG7", "PSMB10"),
  Checkpoint = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA", "ENTPD1", "VSIR"),
  Exhaustion = c("TOX", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "CXCL13", "LAYN", "ENTPD1"),
  MHC_I = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"),
  MHC_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CIITA", "CD74"),
  Myeloid_suppression = c("CD274", "IL10", "TGFB1", "ARG1", "IDO1", "S100A8", "S100A9", "MRC1", "MSR1", "MARCO"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8", "ENTPD1"),
  Proliferation = c("MKI67", "TOP2A", "PCNA", "STMN1", "TYMS", "MCM2", "MCM5"),
  Angiogenesis = c("VEGFA", "KDR", "FLT1", "ANGPT2", "PECAM1", "VWF", "ESM1"),
  TGFb = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "COL1A1", "COL1A2"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22", "BLK"),
  Plasma = c("MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC", "DERL3"),
  NK = c("KLRD1", "NCR1", "NKG7", "GNLY", "PRF1", "GZMB", "FCGR3A"),
  Monocyte_Macrophage = c("LYZ", "CST3", "LST1", "FCGR3A", "MS4A7", "CD68", "CD163", "MRC1"),
  Dendritic = c("ITGAX", "CLEC9A", "XCR1", "FCER1A", "CD1C", "LILRA4", "CLEC10A")
)

single_gene_features <- c(
  "CD8A", "CD8B", "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
  "IFNG", "CXCL9", "CXCL10", "GZMB", "PRF1", "FOXP3", "IDO1", "B2M",
  "HLA-DRA", "VEGFA", "TGFB1", "MKI67"
)

compute_signature_features <- function(expr_mat) {
  log_expr <- log2(expr_mat + 1)
  gene_z <- t(scale(t(log_expr)))
  gene_z[!is.finite(gene_z)] <- 0

  sig_df <- data.frame(Run = colnames(expr_mat), check.names = FALSE)
  for (nm in names(signature_sets)) {
    genes <- intersect(toupper(signature_sets[[nm]]), rownames(gene_z))
    sig_df[[paste0("sig_", nm)]] <- if (length(genes)) colMeans(gene_z[genes, , drop = FALSE]) else NA_real_
    sig_df[[paste0("sig_", nm, "_raw")]] <- if (length(genes)) colMeans(log_expr[genes, , drop = FALSE]) else NA_real_
  }
  for (g in intersect(toupper(single_gene_features), rownames(log_expr))) {
    sig_df[[paste0("gene_", g)]] <- as.numeric(log_expr[g, ])
  }
  sig_mat <- as.matrix(sig_df[, -1, drop = FALSE])
  rownames(sig_mat) <- sig_df$Run
  sig_mat[!is.finite(sig_mat)] <- 0
  sig_mat
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

remove_bad_cols <- function(x_train, x_test) {
  s <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  list(train = x_train[, keep, drop = FALSE], test = x_test[, keep, drop = FALSE])
}

predict_glmnet <- function(y_train, x_train, x_test, alpha = 0) {
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(x_train, y_train, family = "binomial", alpha = alpha,
                               type.measure = "auc", nfolds = 5)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

predict_ranger <- function(y_train, x_train, x_test, num.trees = 800, mtry_frac = 0.35, min.node.size = 8) {
  dat <- data.frame(y = factor(y_train, levels = c(0, 1)), as.data.frame(x_train), check.names = FALSE)
  fit <- tryCatch(
    ranger(
      y ~ ., data = dat, probability = TRUE, num.trees = num.trees,
      mtry = max(1L, floor(ncol(x_train) * mtry_frac)),
      min.node.size = min.node.size, seed = sample.int(1e6, 1)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, data = as.data.frame(x_test))$predictions[, "1"])
}

predict_lightgbm <- function(y_train, x_train, x_test, num_leaves = 7, learning_rate = 0.03, nrounds = 180) {
  dtrain <- lgb.Dataset(data = x_train, label = y_train)
  params <- list(
    objective = "binary",
    metric = "auc",
    verbosity = -1,
    num_threads = min(8L, parallel::detectCores(logical = TRUE)),
    learning_rate = learning_rate,
    num_leaves = num_leaves,
    min_data_in_leaf = 25,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    lambda_l1 = 0.05,
    lambda_l2 = 1
  )
  fit <- tryCatch(
    lgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = -1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(predict(fit, x_test))
}

fit_predict <- function(model, y_train, x_train, x_test) {
  if (model == "glmnet_ridge") return(predict_glmnet(y_train, x_train, x_test, alpha = 0))
  if (model == "glmnet_elastic") return(predict_glmnet(y_train, x_train, x_test, alpha = 0.5))
  if (model == "glmnet_lasso") return(predict_glmnet(y_train, x_train, x_test, alpha = 1))
  if (model == "ranger") return(predict_ranger(y_train, x_train, x_test, mtry_frac = 0.35, min.node.size = 8))
  if (model == "ranger_smallnode") return(predict_ranger(y_train, x_train, x_test, mtry_frac = 0.5, min.node.size = 3))
  if (model == "lightgbm_shallow") return(predict_lightgbm(y_train, x_train, x_test, num_leaves = 7, learning_rate = 0.03, nrounds = 180))
  if (model == "lightgbm_mid") return(predict_lightgbm(y_train, x_train, x_test, num_leaves = 15, learning_rate = 0.02, nrounds = 260))
  stop("Unknown model: ", model)
}

run_repeated_cv <- function(x, y, model, label, repeats = 10, k = 5) {
  auc_rows <- list()
  pred_rows <- list()
  for (rr in seq_len(repeats)) {
    folds <- make_stratified_folds(y, k = k)
    pred <- rep(NA_real_, length(y))
    for (ff in seq_len(k)) {
      tr <- which(folds != ff)
      te <- which(folds == ff)
      xb <- remove_bad_cols(x[tr, , drop = FALSE], x[te, , drop = FALSE])
      xs <- scale_by_train(xb$train, xb$test)
      pred[te] <- fit_predict(model, y[tr], xs$train, xs$test)
    }
    auc_rows[[rr]] <- data.frame(model = label, repeat_id = rr, AUC = auc_manual(y, pred),
                                 Mean_R = mean(pred[y == 1], na.rm = TRUE),
                                 Mean_NR = mean(pred[y == 0], na.rm = TRUE))
    pred_rows[[rr]] <- data.frame(model = label, repeat_id = rr, Run = names(y), y = y,
                                  prediction = pred, stringsAsFactors = FALSE)
    message(label, " repeat ", rr, "/", repeats, " AUC=", round(auc_rows[[rr]]$AUC, 4))
  }
  list(auc = do.call(rbind, auc_rows), pred = do.call(rbind, pred_rows))
}

message("Reading clinical scores and ImmuCellAI2 fractions...")
score_df <- fread(score_file, data.table = FALSE, check.names = FALSE)
frac_df <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
names(frac_df)[1] <- "Run"
rownames(frac_df) <- as.character(frac_df$Run)
frac_df$Run <- NULL

y_all <- clean_response(score_df)
score_df$ResponseBinaryClean <- y_all
score_df$Run <- as.character(score_df$Run)
score_df <- score_df[!is.na(score_df$ResponseBinaryClean) & score_df$Run %in% rownames(frac_df), , drop = FALSE]
frac_mat <- as.matrix(frac_df[score_df$Run, , drop = FALSE])
storage.mode(frac_mat) <- "double"
frac_mat[!is.finite(frac_mat)] <- 0

y <- as.integer(score_df$ResponseBinaryClean)
names(y) <- score_df$Run
message("Samples: ", length(y), "; R=", sum(y == 1), "; NR=", sum(y == 0))

message("Reading expression matrix for matched samples...")
expr_mat <- read_matrix_gene_first(expr_file, select_samples = score_df$Run)
expr_mat <- expr_mat[, score_df$Run, drop = FALSE]
message("Expression genes: ", nrow(expr_mat), "; expression samples: ", ncol(expr_mat))

sig_mat <- compute_signature_features(expr_mat)
sig_mat <- sig_mat[score_df$Run, , drop = FALSE]
fwrite(data.frame(Run = rownames(sig_mat), sig_mat, check.names = FALSE),
       file.path(out_dir, "immunotherapy_expression_signature_features.txt"), sep = "\t")

score_cols <- intersect(c(
  "C3_like_spearman", "C4_like_spearman", "C3_like_pearson", "C4_like_pearson",
  "C3_like_cosine", "C4_like_cosine", "Weighted_C4_minus_C3",
  "C4_minus_C3_spearman", "C4_minus_C3_pearson", "C4_minus_C3_cosine"
), colnames(score_df))
score_mat <- as.matrix(score_df[, score_cols, drop = FALSE])
storage.mode(score_mat) <- "double"
score_mat[!is.finite(score_mat)] <- 0
rownames(score_mat) <- score_df$Run
colnames(score_mat) <- paste0("score_", make.names(colnames(score_mat)))

cell_asin <- asin(sqrt(pmax(frac_mat, 0)))
rownames(cell_asin) <- score_df$Run
colnames(cell_asin) <- paste0("cell_", make.names(colnames(cell_asin)))

cell_logit <- log((pmax(frac_mat, 0) + 1e-5) / (1 - pmin(frac_mat, 1 - 1e-5)))
rownames(cell_logit) <- score_df$Run
colnames(cell_logit) <- paste0("cell_logit_", make.names(colnames(cell_logit)))

meta_cols <- intersect(c("Cancer", "Cancer_type", "Tissue", "Drug", "Anti_target", "Biopsy_Time", "Gender", "disease"), colnames(score_df))
meta_mat <- one_hot_meta(score_df, meta_cols)
rownames(meta_mat) <- score_df$Run

study_cols <- intersect(c(meta_cols, "SRA_study", "SRA_Study", "batch"), colnames(score_df))
study_meta_mat <- one_hot_meta(score_df, study_cols)
rownames(study_meta_mat) <- score_df$Run

feature_list <- list(
  strict_cells_scores = cbind(cell_asin, score_mat),
  strict_cells_scores_signatures = cbind(cell_asin, score_mat, sig_mat),
  clinical_cells_scores = cbind(cell_asin, score_mat, meta_mat),
  clinical_cells_scores_signatures = cbind(cell_asin, score_mat, sig_mat, meta_mat),
  upper_study_cells_scores_signatures = cbind(cell_asin, score_mat, sig_mat, study_meta_mat)
)

grid <- expand.grid(
  feature_set = names(feature_list),
  model = c("glmnet_ridge", "glmnet_elastic", "glmnet_lasso", "ranger", "ranger_smallnode", "lightgbm_shallow", "lightgbm_mid"),
  stringsAsFactors = FALSE
)

all_auc <- list()
all_pred <- list()

for (ii in seq_len(nrow(grid))) {
  fs <- grid$feature_set[ii]
  model <- grid$model[ii]
  label <- paste(fs, model, sep = "__")
  message("Running ", label)
  x <- as.matrix(feature_list[[fs]])
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- 0
  reps <- if (grepl("^upper_study", fs)) 5 else 10
  res <- run_repeated_cv(x, y, model, label, repeats = reps, k = 5)
  all_auc[[label]] <- res$auc
  all_pred[[label]] <- res$pred
}

auc_by_repeat <- do.call(rbind, all_auc)
pred_all <- do.call(rbind, all_pred)

summary_df <- aggregate(cbind(AUC, Mean_R, Mean_NR) ~ model, data = auc_by_repeat, FUN = mean)
summary_df$MedianAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = median)$AUC
summary_df$SDAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = sd)$AUC
summary_df$MaxAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = max)$AUC
summary_df$MinAUC <- aggregate(AUC ~ model, data = auc_by_repeat, FUN = min)$AUC
summary_df <- summary_df[order(summary_df$AUC, decreasing = TRUE), ]

fwrite(auc_by_repeat, file.path(out_dir, "signature_lightgbm_ranger_auc_by_repeat.txt"), sep = "\t")
fwrite(pred_all, file.path(out_dir, "signature_lightgbm_ranger_predictions.txt"), sep = "\t")
fwrite(summary_df, file.path(out_dir, "signature_lightgbm_ranger_auc_summary.txt"), sep = "\t")

best_model <- summary_df$model[1]
best_pred <- pred_all[pred_all$model == best_model & pred_all$repeat_id == 1, , drop = FALSE]

p_auc <- ggplot(auc_by_repeat, aes(x = reorder(model, AUC, median), y = AUC)) +
  geom_boxplot(fill = "#6BA5C3", color = "black", outlier.size = 0.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#F66463") +
  coord_flip() +
  theme_bw(base_size = 8) +
  labs(x = NULL, y = "Repeated 5-fold CV AUC",
       title = "AUC optimization with ImmuCellAI2 fractions and expression signatures") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "signature_lightgbm_ranger_auc_boxplot.pdf"), p_auc, width = 10, height = 10)
ggsave(file.path(out_dir, "signature_lightgbm_ranger_auc_boxplot.png"), p_auc, width = 10, height = 10, dpi = 300)

roc_obj <- tryCatch(pROC::roc(best_pred$y, best_pred$prediction, quiet = TRUE), error = function(e) NULL)
if (!is.null(roc_obj)) {
  pdf(file.path(out_dir, "best_signature_model_ROC_repeat1.pdf"), width = 4.5, height = 4.5)
  plot(roc_obj, col = "#F66463", lwd = 2, legacy.axes = TRUE,
       main = paste0("Best model ROC: AUC=", round(as.numeric(pROC::auc(roc_obj)), 3)))
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  dev.off()
}

plot_df <- best_pred
plot_df$response <- factor(ifelse(plot_df$y == 1, "Responder", "Non-responder"),
                           levels = c("Non-responder", "Responder"))
p_resp <- ggplot(plot_df, aes(x = response, y = prediction, fill = response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
  scale_fill_manual(values = c("Non-responder" = "#95A8AC", "Responder" = "#F66463")) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Predicted response probability",
       title = best_model,
       subtitle = paste0("Repeat 1 AUC = ", round(auc_manual(plot_df$y, plot_df$prediction), 3))) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "best_signature_model_response_boxplot_repeat1.pdf"), p_resp, width = 4.5, height = 4.2)
ggsave(file.path(out_dir, "best_signature_model_response_boxplot_repeat1.png"), p_resp, width = 4.5, height = 4.2, dpi = 300)

writeLines(c(
  paste0("Best model: ", best_model),
  paste0("Mean repeated 5-fold CV AUC: ", signif(summary_df$AUC[1], 5)),
  paste0("Median repeated 5-fold CV AUC: ", signif(summary_df$MedianAUC[1], 5)),
  paste0("Max repeat AUC: ", signif(summary_df$MaxAUC[1], 5)),
  "",
  "Feature sets:",
  "strict_cells_scores = ImmuCellAI2 fractions + C3/C4 scores only.",
  "strict_cells_scores_signatures = strict features + expression immune signatures.",
  "clinical_cells_scores_signatures = strict + signatures + clinical context.",
  "upper_study_cells_scores_signatures = additionally includes study/batch labels; treat as optimistic upper bound."
), con = file.path(out_dir, "signature_lightgbm_ranger_readme.txt"))

print(summary_df)
