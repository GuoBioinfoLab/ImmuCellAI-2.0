suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript analysis/icb/immuicbscore_pooled_fivefold_cv.R",
      "fraction_file clinical_file output_dir [n_threads]"
    )
  )
}

fraction_file <- args[1]
clinical_file <- args[2]
output_dir <- args[3]
n_threads <- if (length(args) >= 4L) as.integer(args[4]) else 1L
if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

selected_studies <- c(
  "anti-PD1_SRP070710",
  "anti-PD1_SRP230414",
  "anti-PD1_SRP351936",
  "anti-PD1_ERP105482",
  "anti-PD1-anti-CTLA4_ERP105482",
  "anti-PD1_ERP107734",
  "anti-PD1_ERP117672",
  "anti-CTLA4-to-anti-PD1_SRP417444"
)

clr_transform <- function(x, eps = 1e-5) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx), "-")
}

fit_scaler <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  scale_value <- apply(x, 2, stats::sd, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  list(center = center, scale = scale_value)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

auc_manual <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2L) return(NA_real_)
  n_positive <- sum(y == 1L)
  n_negative <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

roc_coordinates <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) < 2L) return(NULL)
  order_index <- order(score, decreasing = TRUE)
  y <- y[order_index]
  score <- score[order_index]
  n_positive <- sum(y == 1L)
  n_negative <- sum(y == 0L)
  data.table(
    threshold = c(Inf, score, -Inf),
    FPR = c(0, cumsum(y == 0L) / n_negative, 1),
    TPR = c(0, cumsum(y == 1L) / n_positive, 1)
  )
}

make_stratified_folds <- function(y, k = 5L, seed = 20260701L) {
  set.seed(seed)
  fold <- integer(length(y))
  for (class_value in sort(unique(y))) {
    index <- sample(which(y == class_value))
    fold[index] <- rep(seq_len(k), length.out = length(index))
  }
  fold
}

fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
clinical <- fread(clinical_file, data.table = FALSE, check.names = FALSE)
if (ncol(fraction) != 54L) {
  warning(
    "Expected one sample column plus 53 fractions, but found ",
    ncol(fraction), " columns."
  )
}

names(fraction)[1] <- "sample"
if (!all(c("Run", "SRA_study", "ResponseBinary") %in% names(clinical))) {
  stop("Clinical file must contain Run, SRA_study, and ResponseBinary columns.")
}

clinical_keep <- clinical[, c("Run", "SRA_study", "ResponseBinary"), drop = FALSE]
names(clinical_keep)[1] <- "sample"
dat <- merge(fraction, clinical_keep, by = "sample")
dat <- dat[
  dat$SRA_study %in% selected_studies &
    dat$ResponseBinary %in% c(0, 1) &
    !is.na(dat$ResponseBinary),
  ,
  drop = FALSE
]

cell_columns <- setdiff(names(fraction), "sample")
dat <- dat[complete.cases(dat[, cell_columns, drop = FALSE]), , drop = FALSE]
dat$ResponseBinary <- as.integer(dat$ResponseBinary)

missing_studies <- setdiff(selected_studies, unique(dat$SRA_study))
if (length(missing_studies) > 0L) {
  stop("Selected studies missing after alignment: ", paste(missing_studies, collapse = ", "))
}
if (nrow(dat) != 334L) {
  warning("The manuscript analysis contains 334 samples; aligned input contains ", nrow(dat), ".")
}

folds <- make_stratified_folds(dat$ResponseBinary, k = 5L, seed = 20260701L)
prediction_rows <- vector("list", 5L)

for (fold_id in seq_len(5L)) {
  train <- dat[folds != fold_id, , drop = FALSE]
  held_out <- dat[folds == fold_id, , drop = FALSE]

  train_clr <- clr_transform(train[, cell_columns, drop = FALSE])
  held_out_clr <- clr_transform(held_out[, cell_columns, drop = FALSE])
  scaler <- fit_scaler(train_clr)
  train_x <- apply_scaler(train_clr, scaler)
  held_out_x <- apply_scaler(held_out_clr, scaler)

  train_y <- train$ResponseBinary
  class_weights <- c(
    "0" = 0.5 / mean(train_y == 0L),
    "1" = 0.5 / mean(train_y == 1L)
  )

  forest_seed <- 9001L + fold_id
  set.seed(forest_seed)
  model <- ranger(
    x = as.data.frame(train_x, check.names = FALSE),
    y = factor(train_y, levels = c(0, 1)),
    probability = TRUE,
    classification = TRUE,
    num.trees = 1200,
    mtry = max(1L, floor(sqrt(ncol(train_x)))),
    min.node.size = 8,
    class.weights = class_weights,
    seed = forest_seed,
    num.threads = n_threads,
    importance = "impurity"
  )

  score <- predict(
    model,
    data = as.data.frame(held_out_x, check.names = FALSE),
    num.threads = n_threads
  )$predictions[, "1"]

  prediction_rows[[fold_id]] <- data.frame(
    sample = held_out$sample,
    SRA_study = held_out$SRA_study,
    ResponseBinary = held_out$ResponseBinary,
    fold = fold_id,
    ImmuICBscore = as.numeric(score),
    stringsAsFactors = FALSE
  )
}

predictions <- rbindlist(prediction_rows, fill = TRUE)
setorder(predictions, sample)
pooled_auc <- auc_manual(predictions$ResponseBinary, predictions$ImmuICBscore)

summary_table <- data.table(
  Design = "Eight selected ICB cohorts; pooled stratified five-fold cross-validation",
  N_studies = uniqueN(predictions$SRA_study),
  N = nrow(predictions),
  N_R = sum(predictions$ResponseBinary == 1L),
  N_NR = sum(predictions$ResponseBinary == 0L),
  AUC = pooled_auc
)

per_study <- predictions[, .(
  N = .N,
  N_R = sum(ResponseBinary == 1L),
  N_NR = sum(ResponseBinary == 0L),
  AUC = auc_manual(ResponseBinary, ImmuICBscore),
  Mean_R = mean(ImmuICBscore[ResponseBinary == 1L]),
  Mean_NR = mean(ImmuICBscore[ResponseBinary == 0L])
), by = SRA_study][order(-AUC)]

roc <- roc_coordinates(predictions$ResponseBinary, predictions$ImmuICBscore)
fwrite(predictions, file.path(output_dir, "pooled_fivefold_predictions.tsv"), sep = "\t")
fwrite(summary_table, file.path(output_dir, "pooled_fivefold_summary.tsv"), sep = "\t")
fwrite(per_study, file.path(output_dir, "per_study_auc.tsv"), sep = "\t")
fwrite(roc, file.path(output_dir, "pooled_roc_coordinates.tsv"), sep = "\t")

label <- sprintf(
  "AUC = %.3f\nN = %d\nR = %d, NR = %d",
  pooled_auc,
  nrow(predictions),
  sum(predictions$ResponseBinary == 1L),
  sum(predictions$ResponseBinary == 0L)
)

plot_roc <- ggplot(roc, aes(FPR, TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_path(color = "#B2182B", linewidth = 1) +
  annotate("text", x = 0.53, y = 0.18, label = label, hjust = 0, fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "ImmuICBscore ROC",
    subtitle = "Eight selected cohorts; pooled five-fold cross-validation",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(output_dir, "ImmuICBscore_pooled_ROC.pdf"), plot_roc, width = 4.9, height = 4.6)
ggsave(file.path(output_dir, "ImmuICBscore_pooled_ROC.png"), plot_roc, width = 4.9, height = 4.6, dpi = 400)

print(summary_table)
print(per_study)
