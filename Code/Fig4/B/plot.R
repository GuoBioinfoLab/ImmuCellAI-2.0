# Sourced by ../plot_panels.R after statistical verification.
# Calendar-average conversion is explicit; the archived endpoint times remain days.
kmp <- rbind(data.frame(phenotype=names(ip_col),time_days=0,survival=1),km[,c("phenotype","time_days","survival")])
kmp <- kmp[order(kmp$phenotype,kmp$time_days),]
save_plot(ggplot(kmp,aes(time_days/30.4375,survival,colour=phenotype)) + geom_step(linewidth=.6) +
  scale_colour_manual(values=ip_col) + labs(x="Months (days / 30.4375)",y="Overall survival",colour=NULL) +
  annotate("text",x=230,y=.8,label="P < 0.0001") + base_theme + theme(legend.position="top"),"B",5,4.3)
