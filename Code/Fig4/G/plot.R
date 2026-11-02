# Sourced by ../plot_panels.R after statistical verification.
gd <- transform(gdata,phenotype=sub("C","IP",cluster))
gd$score_quartile <- factor(gd$score_quartile,levels=levels(quartiles))
save_plot(ggplot(gd,aes(score_quartile,Proportion,fill=phenotype)) + geom_col(colour="black",linewidth=.4) +
  scale_fill_manual(values=ip_col) + scale_y_continuous(labels=scales::percent) +
  labs(x="ImmuICBscore quartile",y="Proportion") + base_theme + theme(axis.text.x=element_text(angle=25,hjust=1)),"G",5,4)
