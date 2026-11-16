
tab <- function(p) read_tsv(file.path(data_dir, p))
d <- tab("I/ICB_scores.tsv"); predictions <- tab("I/ICB_repeated_OOF_predictions.tsv")
folds <- tab("I/ICB_repeated_folds.tsv")
stopifnot(nrow(d) == 312L, !anyDuplicated(d$sample), length(unique(d$grouping_id)) == 273L,
          !anyDuplicated(table_key(d, c("grouping_id", "timing"))), nrow(predictions) == 6240L,
          !anyDuplicated(table_key(predictions, c("sample", "repeat_id"))),
          all(table(predictions$sample) == 20L),
          all(tapply(folds$fold, table_key(folds, c("repeat_id", "grouping_id")), function(x) length(unique(x))) == 1))
means <- aggregate(score ~ sample, predictions, mean)
check_num(d$score, means$score[match(d$sample, means$sample)], "Mean OOF score")
result <- tab("I/AUC_all_subgroups.tsv"); draws <- tab("I/bootstrap_draws.tsv.gz"); repeats <- tab("I/repeat_AUC.tsv")
for (i in seq_len(nrow(result))) {
  row <- result[i, ]; sub <- d[d$timing == row$timing, ]
  if (row$dimension != "overall") sub <- sub[sub[[row$dimension]] == row$subgroup, ]
  stopifnot(nrow(sub) == row$n, sum(sub$ResponseBinary) == row$responders)
  if (row$status != "estimated") {
    stopifnot(length(unique(sub$ResponseBinary)) == 1L, is.na(row$AUC)); next
  }
  check_num(auc_rank(sub$ResponseBinary, sub$score), row$AUC)
  z <- draws[draws$timing == row$timing & draws$dimension == row$dimension & draws$subgroup == row$subgroup, ]
  stopifnot(nrow(z) == 2000L, !anyDuplicated(z$bootstrap))
  check_num(quantile(z$AUC, c(.025, .975), type = 7), c(row$CI_low, row$CI_high), "Frozen bootstrap CI")
  for (rep_id in 1:20) {
    p <- predictions[predictions$sample %in% sub$sample & predictions$repeat_id == rep_id, ]
    rr <- repeats[repeats$timing == row$timing & repeats$dimension == row$dimension &
                   repeats$subgroup == row$subgroup & repeats$repeat_id == rep_id, ]
    stopifnot(nrow(rr) == 1L)
    check_num(auc_rank(sub$ResponseBinary[match(p$sample, sub$sample)], p$score), rr$AUC)
  }
}
display <- tab("I/plot_data.tsv")
keys <- c("timing", "dimension", "subgroup")
expected <- result[match(table_key(display, keys), table_key(result, keys)), ]
check_table(display, expected, keys, names(expected))
stopifnot(identical(as.integer(display$n), c(243L,163L,45L,35L,188L,32L,23L,69L)))
analysis$I <- display
