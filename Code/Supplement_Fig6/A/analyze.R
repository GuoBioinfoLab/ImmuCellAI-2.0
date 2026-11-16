
tab <- function(p) read_tsv(file.path(data_dir, p))
summary <- tab("shared/feature_stability/all53_feature_summary.tsv")
imp <- tab("A/repeat_importance_all53.tsv")
ranks <- tab("B/ranks_all53.tsv")
perm <- tab("A/fold_permutation_auc_decreases.tsv")
folds <- tab("shared/feature_stability/repeated_folds.tsv")
pred <- tab("shared/feature_stability/out_of_fold_predictions.tsv")
perf <- tab("shared/feature_stability/fold_performance.tsv")
audit <- tab("shared/feature_stability/sample_patient_audit.tsv")
stopifnot(nrow(summary) == 53L, nrow(imp) == 53L, nrow(perm) == 106000L, nrow(audit) == 334L,
          !anyDuplicated(audit$sample), sum(audit$ResponseBinary) == 120L, nrow(pred) == 6680L, nrow(perf) == 100L,
          !anyDuplicated(table_key(perm, c("repeat_id", "fold", "feature", "permutation"))),
          !anyDuplicated(table_key(folds, c("repeat", "sample"))),
          all(tapply(folds$fold, table_key(folds, c("repeat", "grouping_id")), function(x) length(unique(x))) == 1))
for (i in seq_len(nrow(perf))) {
  z <- perf[i, ]; g <- pred[pred$repeat_id == z$repeat_id & pred$fold == z$fold, ]
  expected <- folds[folds[["repeat"]] == z$repeat_id & folds$fold == z$fold, ]
  stopifnot(setequal(g$sample, expected$sample), nrow(g) == z$n_validation)
  check_num(g$response, expected$ResponseBinary[match(g$sample, expected$sample)])
  point <- auc_rank(g$response, g$probability)
  check_num(point, z$baseline_auc)
  pp <- perm[perm$repeat_id == z$repeat_id & perm$fold == z$fold, ]
  check_num(pp$baseline_auc, rep(point, nrow(pp)))
}
for (rep_id in 1:20) {
  p <- perm[perm$repeat_id == rep_id, ]
  a <- aggregate(auc_decrease ~ feature + fold, p, mean)
  w <- perf[perf$repeat_id == rep_id, ]
  a$weighted <- a$auc_decrease * w$n_validation[match(a$fold, w$fold)] / 334
  a <- aggregate(weighted ~ feature, a, sum)
  v <- a$weighted[match(imp$feature, a$feature)]
  check_num(v, imp[[paste0("repeat", rep_id)]], "Repeated importance")
  check_num(rank(-v, ties.method = "first"), ranks[[paste0("repeat", rep_id)]][match(imp$feature, ranks$feature)], "Ranks")
}
s <- summary[match(imp$feature, summary$feature), ]
values <- as.matrix(imp[-1]); rk <- as.matrix(ranks[match(imp$feature, ranks$feature), -1])
check_num(s$mean_auc_decrease, rowMeans(values)); check_num(s$sd_auc_decrease, apply(values, 1, sd))
check_num(s$top10_count, rowSums(rk <= 10)); check_num(s$top10_frequency, s$top10_count / 20)
stopifnot(sum(s$top10_count) == 200L, all(diff(summary$mean_auc_decrease) <= 0))
top <- head(summary, 15)
stopifnot(identical(top$feature, c("Tex","MAIT","pDC","ILC2","CD4Tem","gdT","CD8Temra","M2","Treg","Bex","M0","CD8Tn","MZB","CD4Tn","M1")))
a <- tab("A/plot_data.tsv"); check_table(a, top, "feature", names(a))
analysis$A <- a
