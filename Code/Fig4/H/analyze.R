# Sourced by ../reproduce_Fig4.R in panel order.
# H: AUC is calculated from the archived out-of-fold predictions, not from a
# refitted forest or the historical cohort-search procedure.
hdata <- read_tab("H/predictions.tsv")
hsummary <- read_tab("H/AUC_summary.tsv")
hroc <- read_tab("H/ROC_plot_data.tsv")
hscores <- read_tab("H/ROC_scores.tsv")
hmetrics <- read_tab("H/model_comparison.tsv")
check("H_sample_counts", nrow(hdata)==334 && sum(hdata$ResponseBinary==1)==120 &&
      sum(hdata$ResponseBinary==0)==214 && !anyDuplicated(hdata$sample))
check("CH_same_samples", setequal(cdata$sample, hdata$sample))
check("H_five_folds", setequal(hdata$fold,1:5))
y <- hdata$ResponseBinary
score <- hdata$Unified_ICB_Score
auc <- (sum(rank(score)[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))
check("H_AUC_0806", near(auc, hsummary$AUC) && near(auc, 0.805607476635514))
hmodels <- c("Fine53_RF", "Canonical12_RF", "CD45_direct")
check("H_three_models", identical(hmetrics$model, hmodels) && setequal(hroc$model, hmodels))
check("H_three_curve_samples", nrow(hscores) == 334 && !anyDuplicated(hscores$sample) &&
      setequal(hscores$sample, hdata$sample))
ix <- match(hscores$sample, hdata$sample)
check("H_labels_and_fine_scores", identical(hscores$ResponseBinary, hdata$ResponseBinary[ix]) &&
      near(hscores$Fine53_RF, hdata$Unified_ICB_Score[ix]))
cd45 <- read_tab("H/CD45_expression.tsv")
check("H_direct_CD45_no_RF", near(hscores$CD45_direct, cd45$PTPRC_TPM[match(hscores$sample, cd45$sample)]))
fractions <- read_tab("H/fine_fractions.tsv")
names(fractions)[1] <- "sample"
fine <- fractions[match(hscores$sample, fractions$sample), , drop=FALSE]
check("H_fine_inputs_334_by_53", nrow(fine) == 334 && ncol(fine) == 54 &&
      identical(fine$sample, hscores$sample))
canonical <- read_tab("H/canonical_fractions.tsv")
mapping <- read_tab("H/lineage_mapping.tsv")
check("H_53_states_mapped_once", nrow(mapping) == 53 && !anyDuplicated(mapping$state) &&
      setequal(mapping$state, names(fine)[-1]) && length(unique(mapping$group)) == 12)
for (g in unique(mapping$group)) {
  expected <- rowSums(fine[, mapping$state[mapping$group == g], drop=FALSE])
  check(paste0("H_group_", g), near(expected, canonical[[g]][match(fine$sample, canonical$sample)]))
}
rocs <- lapply(hmodels, function(m) pROC::roc(hscores$ResponseBinary, hscores[[m]],
                                           levels=c(0,1), direction="<", quiet=TRUE))
names(rocs) <- hmodels
hs <- hmetrics
recomputed_rocs <- do.call(rbind, lapply(hmodels, function(m) {
  rr <- rocs[[m]]
  tab <- data.frame(model=m, FPR=1-rr$specificities, TPR=rr$sensitivities, threshold=rr$thresholds)
  saved <- hroc[hroc$model == m, ]
  check(paste0("H_ROC_", m), near(tab$FPR, saved$FPR) && near(tab$TPR, saved$TPR))
  rank_auc <- (sum(rank(hscores[[m]])[hscores$ResponseBinary==1])-120*121/2)/(120*214)
  check(paste0("H_AUC_", m), near(rank_auc, as.numeric(pROC::auc(rr))) &&
        near(rank_auc, hmetrics$AUC[hmetrics$model==m]))
  hs$AUC[hs$model==m] <<- rank_auc
  tab
}))
