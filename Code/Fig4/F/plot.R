# Sourced by ../plot_panels.R after statistical verification.
fd <- transform(fdata,phenotype=sub("C","IP",cluster))
save_plot(ggplot(fd,aes(phenotype,Unified_ICB_score,fill=phenotype)) +
  geom_violin(trim=TRUE,linewidth=.5) + geom_boxplot(width=.16,outlier.shape=NA,linewidth=.5) +
  stat_summary(fun=mean,geom="point",shape=23,fill="white",size=2) + scale_fill_manual(values=ip_col) +
  labs(x="Immune phenotypes",y="ImmuICBscore") + base_theme,"F",4.5,4)
