# Sourced by ../reproduce_Fig4.R in panel order.
# F/G: freeze existing model scores. The historical TCGA score model was
# trained on all eligible ICB training samples, unlike the H out-of-fold model.
fdata <- read_tab("F/score_plot_data.tsv")
fsummary <- read_tab("F/score_summary.tsv")
check("F_4854_scores_unique", nrow(fdata) == 4854 && !anyDuplicated(fdata$sample))
check("AF_assignments_match", identical(sub("C", "IP", fdata$cluster), ann$phenotype[match(fdata$sample, ann$sample)]))
fs <- do.call(rbind, lapply(paste0("C",1:4), function(k) {
  v <- fdata$Unified_ICB_score[fdata$cluster == k]
  data.frame(cluster=k, N=length(v), MeanScore=mean(v), MedianScore=median(v))
}))
check("F_summary_reproduced", near(fs$MeanScore, fsummary$MeanScore) && near(fs$N, fsummary$N) && near(fs$MedianScore, fsummary$MedianScore))
breaks <- quantile(fdata$Unified_ICB_score, seq(0,1,.25))
breaks[1] <- breaks[1]-1e-12; breaks[5] <- breaks[5]+1e-12
quartiles <- cut(fdata$Unified_ICB_score, breaks, include.lowest=TRUE,
                labels=c("Q1 lowest", "Q2", "Q3", "Q4 highest"))
check("FG_quartile_labels", identical(as.character(quartiles), fdata$score_quartile))
