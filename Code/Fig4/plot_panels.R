# Called by reproduce_Fig4.R --plots; uses already verified in-memory data.
suppressPackageStartupMessages({library(ggplot2); library(ComplexHeatmap); library(grid)})
plot_dir <- file.path(out,"panels")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
ip_col <- setNames(c("#E15759", "#4E79A7", "#F28E2B", "#76B7B2"),paste0("IP",1:4))
base_theme <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="grey85",linetype="dashed",linewidth=.35),
  axis.text=element_text(colour="black"), strip.background=element_rect(fill="#F5E6E8"),
  legend.position="none", plot.title=element_text(size=12,face="bold"))
save_plot <- function(p,n,w,h) {
  ggsave(file.path(plot_dir,paste0(n,".png")),p,width=w,height=h,dpi=180,bg="white")
  ggsave(file.path(plot_dir,paste0(n,".pdf")),p,width=w,height=h,
         device=grDevices::pdf,useDingbats=FALSE,bg="white")
}

for (panel in c("A","B","C","D","F","G","H")) {
  source(file.path(code_dir,panel,"plot.R"), local=TRUE)
}
cat("Verified-data panel redraws written to:",plot_dir,"\n")
