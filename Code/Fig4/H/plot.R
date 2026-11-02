# Sourced by ../plot_panels.R after statistical verification.
hlabels <- c(Fine53_RF="ImmuICBscore (fine immune subtypes, RF)",
             Canonical12_RF="Major immune cell types (12, RF)",
             CD45_direct="CD45 expression (no RF)")
hlegend <- setNames(sprintf("%s\nAUC %.3f", hlabels[hmodels], hs$AUC[match(hmodels,hs$model)]), hmodels)
save_plot(ggplot(hroc,aes(FPR,TPR,colour=model)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60") +
  geom_path(linewidth=.7) + coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE) +
  scale_colour_manual(values=c(Fine53_RF="#B13E4B",Canonical12_RF="#267F7A",CD45_direct="#5376A0"),
                      breaks=hmodels,labels=hlegend) +
  annotate("text",x=.53,y=.48,hjust=0,label="N = 334\nR = 120, NR = 214",size=3.5,fontface="bold") +
  labs(title="ImmuICBscore vs immune baselines",x="False positive rate",y="True positive rate",colour=NULL) +
  base_theme + theme(legend.position="bottom",legend.text=element_text(size=10)) +
  guides(colour=guide_legend(ncol=1)),"H",5.5,6.2)
