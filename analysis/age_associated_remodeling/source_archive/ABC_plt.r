load("<LOCAL_LEGACY_PROJECT_ROOT>/Figure3/result_ABC_Figure3.Rdata")
long_data2[1:3,1:3] 
unique(long_data2$cell_type)
long_data2$cell_type <- as.character(long_data2$cell_type) 
long_data2[long_data2['cell_type'] == 'Bnaive', ]$cell_type = "Bn"
long_data2[long_data2['cell_type'] == 'CD4Tnaive', ]$cell_type = "CD4Tn"
long_data2[long_data2['cell_type'] == 'CD8Tnaive', ]$cell_type = "CD8Tn"
long_data2[long_data2['cell_type'] == 'CMonocyte', ]$cell_type = "cMo"
long_data2[long_data2['cell_type'] == 'GC_B', ]$cell_type = "BGC"
long_data2[long_data2['cell_type'] == 'IMonocyte', ]$cell_type = "intMo"
long_data2[long_data2['cell_type'] == 'NMonocyte', ]$cell_type = "ncMo"
long_data2[long_data2['cell_type'] == 'basophils', ]$cell_type = "Basophil"
long_data2[long_data2['cell_type'] == 'cytotoxicNK', ]$cell_type = "cNK"
long_data2[long_data2['cell_type'] == 'eosinophils', ]$cell_type = "Eosinophil"
long_data2[long_data2['cell_type'] == 'exhaustedB', ]$cell_type = "Bex"
long_data2[long_data2['cell_type'] == 'exhausted_T', ]$cell_type = "Tex"
long_data2[long_data2['cell_type'] == 'gdT', ]$cell_type = "γδT"
long_data2[long_data2['cell_type'] == 'mast_cell', ]$cell_type = "Mast cell"
long_data2[long_data2['cell_type'] == 'memoryB', ]$cell_type = "MBC"
long_data2[long_data2['cell_type'] == 'monoDC', ]$cell_type = "moDC"
long_data2[long_data2['cell_type'] == 'neutrophils', ]$cell_type = "Neutrophil"
long_data2[long_data2['cell_type'] == 'plasma', ]$cell_type = "PC"
long_data2[long_data2['cell_type'] == 'plasmablast', ]$cell_type = "PB"
long_data2[long_data2['cell_type'] == 'regulatoryNK', ]$cell_type = "NKreg"
unique(long_data2$group)
# long_data2 <- long_data2 %>%
#   mutate(group = case_when(
#     group == "group0to1" ~ "age0-1",
#     group == "group10+" ~ "age10-20",
#     group == "group20+" ~ "age20-30",
#     group == "group30+50" ~ "age30-50",
#     group == "group50+" ~ "age50-70",
#     group %in% c("group70+","group90+") ~ "age70+",
#     TRUE ~ group  # 保持其他值不变
#   )) 
long_data2 <- long_data2 %>%
  mutate(group = case_when(
    group == "group0to1" ~ "age0-1",
    group == "group10+" ~ "age10-20",
    group == "group20+" ~ "age20-30",
    group == "group30+50" ~ "age30-50",
    group == "group50to70" ~ "age50-70",
    group == "group70+" ~ "age70+",
    TRUE ~ group  # 保持其他值不变
  )) 
CD4Tnaive <- long_data2[which(long_data2$cell_type == "CD4Tn"),] 
head(CD4Tnaive) 
p1 = ggplot(CD4Tnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") + 
  ylim(0, 0.37) +  # 设置y轴范围 
  labs(title = "CD4Tn", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(plot.title = element_text(size = 16),
                          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
                          axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p1 
# + stat_compare_means(
#   comparisons = list(c("group0to1", "group10+"),c("group10+", "group20+"),
#                      c("group50+","group70+")), 
#   method = "wilcox.test",  # 使用t检验
#   label = "p.format", # 显示显著性符号
#   label.y = 0.3,
#   size = 5) 

CD8Tnaive <- long_data2[which(long_data2$cell_type == "CD8Tn"),]
CD8Tnaive[1:5,]
p2 = ggplot(CD8Tnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.09) +  # 设置y轴范围 
  labs(title = "CD8Tn", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(plot.title = element_text(size = 16),
                          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
                          axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) +geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                 size=1,se=FALSE, #添加置信区间             
                 color="black",              
                 linetype="solid",              
                 aes(group=1)) 
p2 
# + stat_compare_means( 
# comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
# method = "wilcox.test",  # 使用t检验
# label = "p.format", # 显示显著性符号
# label.y = 0.1, 
# size = 5) 

CMonocyte <- long_data2[which(long_data2$cell_type == "cMo"),]  
head(CMonocyte) 
p3 = ggplot(CMonocyte, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.5) +  # 设置y轴范围 
  labs(title = "cMo", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x,    
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p3 
# + stat_compare_means(
#   comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
#   method = "wilcox.test",  # 使用t检验
#   label = "p.format", # 显示显著性符号
#   label.y = 0.15,
#   size = 7
# ) 
M0 <- long_data2[which(long_data2$cell_type == "M0"),] 
p4 = ggplot(M0, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.05) +  # 设置y轴范围 
  labs(title = "M0", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, 
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p4 

monoDC <- long_data2[which(long_data2$cell_type == "monoDC"),] 
p9 = ggplot(monoDC, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.08) +  # 设置y轴范围 
  labs(title = "monoDC", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, 
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p9
library(cowplot)
combined_plot1 <- plot_grid(p1, p2, p3, p4, ncol = 4)
combined_plot1 

cytotoxicNK <- long_data2[which(long_data2$cell_type == "cNK"),] 
p5 = ggplot(cytotoxicNK, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.4) +  # 设置y轴范围 
  labs(title = "cNK", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.3,
    size = 7) + geom_smooth(method="loess",formula = y ~ x, 
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 
p5 

Tex <- long_data2[which(long_data2$cell_type == "Tex"),] 
p6 = ggplot(Tex, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.16) +  # 设置y轴范围 
  labs(title = "Tex", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.15,
    size = 7) + geom_smooth(method="loess",formula = y ~ x, 
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 
p6 

IMonocyte <- long_data2[which(long_data2$cell_type == "intMo"),] 
p7=ggplot(IMonocyte, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.07) +  
  labs(title = "intMo", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p7 

M1 <- long_data2[which(long_data2$cell_type == "M1"),] 
p8 = ggplot(M1, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.014) +  # 设置y轴范围 
  labs(title = "M1", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p8 

Tc <- long_data2[which(long_data2$cell_type == "Tc"),] 
p10 = ggplot(Tc, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.22) +  
  labs(title = "Tc", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  
  ) + geom_smooth(method="loess",formula = y ~ x, 
                  size=1,se=FALSE, 
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p10 

combined_plot2 <- plot_grid(p5, p6, p7, p8, ncol = 4)
combined_plot2

basophils <- long_data2[which(long_data2$cell_type == "Basophil"),] 
p11 = ggplot(basophils, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.013) +  # 设置y轴范围 
  labs(title = "Basophil", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.15,
    size = 7) + geom_smooth(method="loess",formula = y ~ x, 
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 
p11

Bnaive <- long_data2[which(long_data2$cell_type == "Bn"),] 
p12 = ggplot(Bnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.021) +  # 设置y轴范围 
  labs(title = "Bn", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.15,
    size = 7) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))
p12

Breg <- long_data2[which(long_data2$cell_type == "Breg"),] 
p13 = ggplot(Breg, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.09) +  # 设置y轴范围 
  labs(title = "Breg", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.15,
    size = 7) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))
p13

CD4Tem <- long_data2[which(long_data2$cell_type == "CD4Tem"),] 
head(CD4Tem) 
p14=ggplot(CD4Tem, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.013) +  # 设置y轴范围 
  labs(title = "CD4Tem", x = "Group", y = "Cell Ratios") + 
  theme_minimal()  + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p14

CD8Tem <- long_data2[which(long_data2$cell_type == "CD8Tem"),] 
head(CD8Tem) 
p15 = ggplot(CD8Tem, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.013) +  # 设置y轴范围 
  labs(title = "CD8Tem", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  
  ) +
  geom_smooth(method="loess",formula = y ~ x, 
              size=1,se=FALSE, #添加置信区间             
              color="black",              
              linetype="solid",              
              aes(group=1)) 
p15 

neutrophils <- long_data2[which(long_data2$cell_type == "neutrophils"),] 
p16 = ggplot(neutrophils, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.6) +  # 设置y轴范围 
  labs(title = "neutrophils", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p16

FOB <- long_data2[which(long_data2$cell_type == "FOB"),] 
p17 = ggplot(FOB, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.061) +  # 设置y轴范围 
  labs(title = "FOB", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  )+ geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                 size=1,se=FALSE, #添加置信区间             
                 color="black",              
                 linetype="solid",              
                 aes(group=1))
p17

BGC <- long_data2[which(long_data2$cell_type == "BGC"),] 
p18 = ggplot(BGC, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.012) +  # 设置y轴范围 
  labs(title = "BGC", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p18

gdT <- long_data2[which(long_data2$cell_type == "γδT"),] 
p19 = ggplot(gdT, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.015) +  # 设置y轴范围 
  labs(title = "γδT", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x,  
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p19 

ILC1 <- long_data2[which(long_data2$cell_type == "ILC1"),] 
head(ILC1) 
p20=ggplot(ILC1, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.01125) +  # 设置y轴范围 
  labs(title = "ILC1", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, 
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p20

MAIT <- long_data2[which(long_data2$cell_type == "MAIT"),] 
p21 = ggplot(MAIT, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.012) +  # 设置y轴范围 
  labs(title = "MAIT", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p21

memoryB <- long_data2[which(long_data2$cell_type == "memoryB"),] 
p22 = ggplot(memoryB, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.019) +  # 设置y轴范围 
  labs(title = "memoryB", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p22

MZB <- long_data2[which(long_data2$cell_type == "MZB"),] 
p23 = ggplot(MZB, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.26) +  # 设置y轴范围 
  labs(title = "MZB", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p23

NKT <- long_data2[which(long_data2$cell_type == "NKT"),] 
p24 = ggplot(NKT, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.09) +  # 设置y轴范围 
  labs(title = "NKT", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p24

Th1 <- long_data2[which(long_data2$cell_type == "Th1"),] 
p25 = ggplot(Th1, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.014) +  # 设置y轴范围 
  labs(title = "Th1", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p25

Th2 <- long_data2[which(long_data2$cell_type == "Th2"),] 
p26 = ggplot(Th2, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.0125) +  # 设置y轴范围 
  labs(title = "Th2", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))
p26

Tfh <- long_data2[which(long_data2$cell_type == "Tfh"),] 
p27 = ggplot(Tfh, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.15) +  # 设置y轴范围 
  labs(title = "Tfh", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, 
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1)) 
p27

pDC <- long_data2[which(long_data2$cell_type == "pDC"),] 
p48 =ggplot(pDC, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.013) +  # 设置y轴范围 
  labs(title = "pDC", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))  

combined_plot3 <- plot_grid(p11, p12, p13, p14, p15,
                            p18, p19, p20, p21) 
# p22, p23, p24, p25,
# p26, p27,p48
combined_plot3

combined_all1 <- plot_grid(combined_plot1, 
                           combined_plot2,ncol = 1)
combined_all2 <- plot_grid(combined_plot3, 
                           pE,nrow = 1, rel_widths = c(2.9, 1))
combined_all <- plot_grid(combined_all1, 
                          combined_all2,ncol = 1) 
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/Figure3ABC.pdf",
       combined_all, width = 15, height = 15, device = cairo_pdf) 
save(combined_plot3,long_data2, file = "<LOCAL_LEGACY_PROJECT_ROOT>/Figure3/result_ABC_Figure3.Rdata")
load("<LOCAL_LEGACY_PROJECT_ROOT>/Figure3/result_ABC_Figure3.Rdata") 

