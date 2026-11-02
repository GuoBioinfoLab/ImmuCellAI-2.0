# Sourced by ../plot_panels.R after statistical verification.
cp <- cdata[retained,]
cells <- c("cDC1","Tc","CD8Tem","MBC")
cp$cell_label <- factor(cp$cell_label,levels=cells)
cp$Response <- factor(cp$Response,levels=c("NR","R"))
cp$Biopsy_Time <- factor(cp$Biopsy_Time,levels=c("Pre","On"))
cp$facet <- factor(paste(cp$cell_label,cp$Response,sep=" | "),levels=c(outer(cells,c("NR","R"),paste,sep=" | ")))
cs <- cstats
cs$facet <- factor(paste(sub("Memory B cell","MBC",cs$Cell_Type),cs$Response,sep=" | "),levels=levels(cp$facet))
cs$height <- vapply(as.character(cs$facet),function(k) max(cp$Cell_Ratio[cp$facet==k])*.94,numeric(1))
dist_layers <- list(geom_violin(trim=TRUE,scale="width",width=.72,alpha=.95,colour=NA),
                    geom_boxplot(width=.42,outlier.shape=NA,alpha=.65,linewidth=.4))
save_plot(ggplot(cp,aes(Biopsy_Time,Cell_Ratio,fill=Biopsy_Time)) + dist_layers +
  scale_fill_manual(values=c(Pre=unname(ip_col[4]),On=unname(ip_col[3]))) + facet_wrap(~facet,ncol=4,scales="free_y") +
  geom_text(data=cs,aes(x=1.5,y=height,label=signif),inherit.aes=FALSE,fontface="bold") +
  labs(x=NULL,y="Cell ratio",title="Pre- and on-treatment differences within response groups") + base_theme,"C",11,6)
