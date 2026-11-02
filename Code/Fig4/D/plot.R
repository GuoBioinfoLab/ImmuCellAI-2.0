# Sourced by ../plot_panels.R after statistical verification.
dp <- cp[cp$Biopsy_Time=="On",]
ds <- dstats
ds$cell_label <- factor(sub("Memory B cell","MBC",ds$Cell_Type),levels=cells)
ds$height <- vapply(as.character(ds$cell_label),function(k) max(dp$Cell_Ratio[dp$cell_label==k])*.94,numeric(1))
save_plot(ggplot(dp,aes(Response,Cell_Ratio,fill=Response)) + dist_layers +
  scale_fill_manual(values=c(NR=unname(ip_col[2]),R=unname(ip_col[1]))) + facet_wrap(~cell_label,nrow=1,scales="free_y") +
  geom_text(data=ds,aes(x=1.5,y=height,label=signif),inherit.aes=FALSE,fontface="bold") +
  labs(x=NULL,y="Cell ratio",title="Response-group differences in on-treatment samples") + base_theme,"D",11,3.2)
