(load("data_long_A.Rdata"))
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggpubr)
library(gghalves) # Compose Half-Half Plots Using Your Favourite Geoms
source("mytheme.R")

data_long$CellType %<>% as.character()

FUN <- function(cell) {
  tmp <- data_long %>% dplyr::filter(CellType == cell)
  ggplot(tmp, aes(x = status, y = Abundance, fill = status)) +
    geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.5, color = "black", size = 0.75) +
    stat_boxplot(geom = "errorbar", width = 0.2, size = 0.75) +
    geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
    stat_summary(fun.y = "median",geom = "point",shape = 16, 
                 size = 2, color = "white",position =position_dodge(1)) +
    # 颜色设置
    scale_fill_manual(values = c("#E15759FF", "#4E79A7FF")) +
    labs(x = NULL, y = "Cell Ratio", title = cell) +
    # 添加显著性标记
    theme(legend.position = "none") +
    theme_blue +
    stat_compare_means(aes(group = status), size = 9, label = "p.signif", hide.ns = F, vjust = 1, label.x = 1.5, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))
}
l1 <- lapply(unique(data_long$CellType) %>% str_sort, FUN)
