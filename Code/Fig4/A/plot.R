# Sourced by ../plot_panels.R after statistical verification.
save_heatmap <- function(format) {
  if (format == "png") {
    png(file.path(plot_dir,"A.png"),width=1500,height=2000,res=180)
  } else {
    pdf(file.path(plot_dir,"A.pdf"),width=1500/180,height=2000/180,
        useDingbats=FALSE,bg="white")
  }
  on.exit(dev.off(), add=TRUE)
ha <- HeatmapAnnotation(phenotype=anno_block(gp=gpar(fill=ip_col), labels=names(ip_col), labels_gp=gpar(col="white")))
heat <- Heatmap(z,name="Row z-score",cluster_rows=FALSE,cluster_columns=FALSE,
                show_column_names=FALSE,column_split=factor(ann$phenotype,levels=names(ip_col)),
                column_title=NULL, top_annotation=ha, col=colorRampPalette(c("#3341A3","black","#E0DA54"))(50),
                row_names_gp=gpar(fontsize=9),border=TRUE,use_raster=TRUE)
draw(heat,heatmap_legend_side="right")
}
save_heatmap("png")
save_heatmap("pdf")
