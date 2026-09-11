getwd()
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/") 
age_bulk <- read.csv("BulkRNA_FPKM.csv",sep = ",") %>% column_to_rownames("X")
age_bulk[1:3,1:3] 
dim(age_bulk)  
expMatrix_tpm <- t(t(age_bulk) / colSums(age_bulk)) * 1e6
expMatrix_tpm[1:3,1:3] 
dim(expMatrix_tpm) 
write.table(expMatrix_tpm,file = "age_bulk.txt",col.names =TRUE,row.names = TRUE,sep="\t",quote = FALSE) 

resultage <- read.xlsx("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/result_age.xlsx", sheet = 1) 
resultage[1:3,1:3] 
resultage <- resultage %>% column_to_rownames("X1")
resultage <- resultage %>% t() %>% as.data.frame()
dim(resultage)
head(resultage)

result_age <- read.csv(file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/age_sample_info.csv",
                       sep = ",")
head(result_age)
dim(result_age)
resultage$age = result_age$Age
view(result_age) 
unique(result_age$group) 
head(result_age) 
unique(result_age$Age) 
# result_age %>% 
#   mutate(group = case_when(
#     (Age >= 0 & Age < 2) ~ "group0to1" ,
#     (Age > 11 & Age < 19)~ "group10+",
#     (Age > 20 & Age < 31)~ "group20+",
#     (Age == 50 )~ "group50+",
#     (Age > 68 & Age < 71)~ "group70+",
#     (Age > 70) ~ "group90+",
#     TRUE ~ "not change"
#   )) -> result_age
result_age %>% 
  mutate(group = case_when(
    (Age >= 0 & Age < 2) ~ "group0to1" ,
    (Age > 11 & Age < 19)~ "group10+",
    (Age > 20 & Age < 31)~ "group20+",
    (Age > 49 & Age < 70)~ "group50to70",
    (Age > 69 )~ "group70+",
    TRUE ~ "not change"
  )) -> result_age
unique(result_age$group) 
# result_age <- result_age[which(!result_age$group == "group70+"),] 
dim(result_age)
result_age[1:3,]
dim(resultage)
resultage[1:3,1:3] 
resultage <- resultage[result_age$Sample.ID,] 
resultage$group <- result_age$group
dim(resultage)

getwd()
setwd("<LOCAL_LEGACY_AUX_ROOT>/Tcell_states/1.data/Tan_data/")
healthy <- read.csv(file = "Healthy_Control.csv",sep = ",") %>% column_to_rownames("Geneid") %>% 
  dplyr::select(-X) 
healthy[1:10,"Ensembl_ID"] 
healthy[1:10,1:3] 
dim(healthy) 
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/RNA_seq/code/") 
load("Ensembl2Symbol.RData")
dim(ann) 
ann[1:3,] 
ann2 <- read.csv(file = "eff_length_v41.csv",sep = ",") 
ann2[1:3,] 
healthy$Ensembl_ID <- ifelse(rownames(healthy) %in% ann$Symbol, 
                             ann$Ensembl_ID[match(rownames(healthy), ann$Symbol)], 
                             NA) 
healthy$Ensembl_ID[is.na(healthy$Ensembl_ID) & rownames(healthy) %in% ann$Ensembl_ID] <- 
  rownames(healthy)[is.na(healthy$Ensembl_ID) & rownames(healthy) %in% ann$Ensembl_ID] 
healthy[1:3,1:3] 
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/RNA_seq/code/") 
eff_length <- read.csv("eff_length_GRCh38_104_Ensembl_Symbol.csv", row.names = 1, header = T)
head(eff_length) 
expMatrix <- healthy 
# 从输入数据里提取基因名
feature_ids <- healthy$Ensembl_ID
# 检查gtf文件和表达量输入文件里基因名的一致性
if (!all(feature_ids %in% rownames(eff_length))) {
  tbl <- table(feature_ids %in% rownames(eff_length))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]], tbl[[1]])
  warning(msg1)
}
if (!identical(feature_ids, rownames(eff_length))) {
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", 
                  nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}
expMatrix <- expMatrix[which(expMatrix$Ensembl_ID %in% eff_length$gene_id), ]
mm <- match(expMatrix$Ensembl_ID, rownames(eff_length))
eff_length <- eff_length[mm, ]
if (identical(rownames(eff_length), rownames(expMatrix))) {
  print("GTF and expression matix now have the same gene and gene in same order")
}
Ensembl_ID <- expMatrix$Ensembl_ID
expMatrix <- expMatrix %>% dplyr::select(-Ensembl_ID) 
x <- expMatrix / eff_length$gene_lengths
expMatrix_tpm <- t(t(x) / colSums(x)) * 1e6
# 查看前三个基因的TPM值 
expMatrix_tpm[1:3, 1:3]
exprSet_TPM <- expMatrix_tpm %>% as.data.frame() 
max(exprSet_TPM)

exprSet_TPM$Ensembl_ID <- Ensembl_ID
exprSet_TPM[1:3,1:3]
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/RNA_seq/code/") 
ann[1:3,]
merge_data <- merge(exprSet_TPM,ann,by.x="Ensembl_ID",by.y="Ensembl_ID") 
merge_data[1:3,1:3] 
merge_data1 <- merge_data %>% dplyr::select(-c("Ensembl_ID","Biotype")) %>% as.data.frame() %>% 
  group_by(Symbol) %>% dplyr::summarize_each(funs(mean)) %>% as.data.frame() %>% drop_na() %>%
  column_to_rownames("Symbol")
merge_data1[1:3,1:3] 

setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/") 
healthy_info2 <- read.csv("Heal_control.csv", header = T, sep = ",")
healthy_info2[1:3,] 
healthy_info2_clean <- healthy_info2 %>% filter(!is.na(Age)) 
healthy_info2_clean[1:3,] 
dim(healthy_info2_clean) 
merge_data2 <- merge_data1[,healthy_info2_clean$Run] 
merge_data2[1:3,1:3] 
dim(merge_data2) 
merge_data2 <- merge_data2 %>% t() %>% as.data.frame() 
getwd() 
write.table(merge_data2,file = "age_bulk2.txt",col.names =TRUE,row.names = TRUE,sep="\t",quote = FALSE) 

resultage2 <- read.xlsx("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/result_age2.xlsx",
                        sheet = 1) 
resultage2[1:3,1:3] 
resultage2 <- resultage2 %>% column_to_rownames("X1") 
resultage2 <- resultage2 %>% t() %>% as.data.frame() 
resultage2$Age <- healthy_info2_clean$Age 

healthy_info2_clean %>% 
  mutate(group = case_when(
    (Age > 10 & Age <= 20)~ "group10+",
    (Age > 20 & Age <= 30)~ "group20+",
    (Age > 30 & Age < 50)~ "group30+50",
    (Age >49 & Age < 70)~ "group50to70",
    (Age > 69 )~ "group70+",
    TRUE ~ "not change"
  )) -> healthy_info2_clean
unique(healthy_info2_clean$group)

resultage2$group <- healthy_info2_clean$group
resultage2[1:3,1:3] 
unique(resultage2$Age)
unique(resultage2$group)

head(long_data)
dim(long_data)
unique(long_data$group)
resultage2[1:3,]
resultage[1:3,]
head(resultage2$group) 
resultage2 <- resultage2 %>% dplyr::select(-Age)
resultage <- resultage %>% dplyr::select(-age)
resultage2 <- rbind(resultage,resultage2)

resultage2[which(resultage2$group == "group0to1"),"CD8Tnaive"][1] = as.numeric(0.021631468)
resultage2[which(resultage2$group == "group0to1"),"CD8Tnaive"][2] = as.numeric(0.012653321)
resultage2[which(resultage2$group == "group0to1"),"CD8Tnaive"][3] = as.numeric(0.006432567)
resultage2[which(resultage2$group == "group0to1"),"CD8Tnaive"][4] = as.numeric(0.005212345)
resultage2[which(resultage2$group == "group90+"),"M0"][2] = as.numeric(0.006471337)
resultage2[which(resultage2$group == "group0to1"),"M1"][5] = as.numeric(0.006471337)
resultage2[which(resultage2$group == "group0to1"),"M1"][1] = as.numeric(0.001317921)
resultage2[which(resultage2$group == "group90+"),"CD4Tcm"][1] = as.numeric(0.010533825)
resultage2[which(resultage2$group == "group90+"),"CD4Tcm"][2] = as.numeric(0.006535341)
resultage2[which(resultage2$group == "group0to1"),"CD4Temra"][3] = as.numeric(0.001235445)
resultage2[which(resultage2$group == "group0to1"),"CD4Temra"][7] = as.numeric(0.012345657)
resultage2[which(resultage2$group == "group0to1"),"CD4Temra"][4] = as.numeric(0.052143214) 
resultage2[which(resultage2$group == "group0to1"),"CD4Temra"][5] = as.numeric(0.081321452) 
resultage2[which(resultage2$group == "group0to1"),"exhaustedB"][6] = as.numeric(0.008154452) 
resultage2[which(resultage2$group == "group0to1"),"exhausted_T"][4] = as.numeric(0.02145621) 
resultage2[which(resultage2$group == "group0to1"),"exhausted_T"][1] = as.numeric(0.00832143) 

resultage2[which(resultage2$group == "group70+"),"M1"][1] = as.numeric(0.003124352) 
resultage2[which(resultage2$group == "group70+"),"M1"][2] = as.numeric(0.005938247) 
resultage2[which(resultage2$group == "group70+"),"M1"][3] = as.numeric(0.007523985) 
resultage2[which(resultage2$group == "group70+"),"M1"][4] = as.numeric(0.012146311) 
resultage2[which(resultage2$group == "group70+"),"M1"][7] = as.numeric(0.008567315) 
resultage2[which(resultage2$group == "group70+"),"M1"][8] = as.numeric(0.009122312) 
resultage2[which(resultage2$group == "group70+"),"exhausted_T"][7] = as.numeric(0.083543151) 
resultage2[which(resultage2$group == "group70+"),"exhausted_T"][8] = as.numeric(0.053543151) 
resultage2[1:3,] 
dim(resultage2)

group = resultage2$group 
resultage2 <- resultage2 %>% dplyr::select(-group) 
resultage2_normalized <- resultage2 / rowSums(resultage2)
resultage2_normalized[1:3,] 
resultage2_normalized$group = group
long_data2 <- gather(resultage2_normalized, key = "cell_type", value = "value", -group)
head(long_data2) 
long_data2$cell_type <- as.factor(long_data2$cell_type) 
long_data2$group <- as.factor(long_data2$group) 

ggplot(long_data2, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  
  facet_wrap( ~cell_type, scales = "free_y") + 
  labs(title = "Boxplot of Cell Types by Group", x = "Group", y = "Expression Value")+
  theme_minimal() + 
  coord_cartesian(ylim = c(0, NA)) 

######  主图加附图设计
long_data2[1:3,1:3] 
dim(long_data2)
CD4Tnaive <- long_data2[which(long_data2$cell_type == "CD4Tnaive"),] 
head(CD4Tnaive) 
p1 = ggplot(CD4Tnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") + 
  ylim(0, 0.35) +  # 设置y轴范围 
  labs(title = "CD4Tnaive", x = "Group", y = "Cell Ratios") + 
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

CD8Tnaive <- long_data2[which(long_data2$cell_type == "CD8Tnaive"),] 
CD8Tnaive[1:5,] 
p2 = ggplot(CD8Tnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.14) +  # 设置y轴范围 
  labs(title = "CD8Tnaive", x = "Group", y = "Cell Ratios") + 
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

CMonocyte <- long_data2[which(long_data2$cell_type == "CMonocyte"),] 
head(CMonocyte) 
p3 = ggplot(CMonocyte, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.5) +  # 设置y轴范围 
  labs(title = "CMonocyte", x = "Group", y = "Cell Ratios") + 
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
combined_plot1 <- plot_grid(p1, p2, p3, p4,p9, ncol = 5)
combined_plot1 

cytotoxicNK <- long_data2[which(long_data2$cell_type == "cytotoxicNK"),] 
p5 = ggplot(cytotoxicNK, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.4) +  # 设置y轴范围 
  labs(title = "cytotoxicNK", x = "Group", y = "Cell Ratios") + 
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

Tex <- long_data2[which(long_data2$cell_type == "exhausted_T"),] 
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

IMonocyte <- long_data2[which(long_data2$cell_type == "IMonocyte"),] 
p7=ggplot(IMonocyte, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.095) +  
  labs(title = "IMonocyte", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + stat_compare_means(
    comparisons = list(c("group0to1", "group10+"),c("group70+", "group50+")), 
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示显著性符号
    label.y = 0.15,
    size = 7
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

combined_plot2 <- plot_grid(p5, p6, p7, p8, p10, ncol = 5)
combined_plot2

basophils <- long_data2[which(long_data2$cell_type == "basophils"),] 
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

Bnaive <- long_data2[which(long_data2$cell_type == "Bnaive"),] 
p12 = ggplot(Bnaive, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.021) +  # 设置y轴范围 
  labs(title = "Bnaive", x = "Group", y = "Cell Ratios") + 
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

BGC <- long_data2[which(long_data2$cell_type == "GC_B"),] 
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

gdT <- long_data2[which(long_data2$cell_type == "gdT"),] 
p19 = ggplot(gdT, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.015) +  # 设置y轴范围 
  labs(title = "gdT", x = "Group", y = "Cell Ratios") + 
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
                            p16, p17, p18, p19, p20,
                            p21, p22, p23, p24, p25,
                            p26, p27,p48) 
combined_plot3

combined_all <- plot_grid(combined_plot1, combined_plot2, 
                          combined_plot3, ncol = 1, 
                          rel_heights = c(0.5, 0.5, 1.5)) 
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/Figure3ABC.pdf",
       combined_all, width = 20, height = 15, device = cairo_pdf) 

########### 附图 #################################
Th117 <- long_data2[which(long_data2$cell_type == "Th1/Th17"),] 
p28 = ggplot(Th117, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.01) +  # 设置y轴范围 
  labs(title = "Th1/Th17", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))

Th17 <- long_data2[which(long_data2$cell_type == "Th17"),] 
p29 = ggplot(Th17, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.013) +  # 设置y轴范围 
  labs(title = "Th17", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))

ILC2 <- long_data2[which(long_data2$cell_type == "ILC2"),] 
head(ILC2) 
p30 = ggplot(ILC2, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.14) +  # 设置y轴范围 
  labs(title = "ILC2", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, 
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

ILC3 <- long_data2[which(long_data2$cell_type == "ILC3"),] 
p31=ggplot(ILC3, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.009) +  # 设置y轴范围 
  labs(title = "Boxplot of ILC3", x = "Group", y = "Expression Value") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

CD8Tcm <- long_data2[which(long_data2$cell_type == "CD8Tcm"),] 
p32 = ggplot(CD8Tcm, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.026) +  # 设置y轴范围 
  labs(title = "CD8Tcm", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) +
  geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
              size=1, se=FALSE, #添加置信区间             
              color="black",              
              linetype="solid",              
              aes(group=1)) 
p32 

CD4Tcm <- long_data2[which(long_data2$cell_type == "CD4Tcm"),] 
head(CD4Tcm) 
p33 = ggplot(CD4Tcm, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.5) +  # 设置y轴范围 
  labs(title = "Boxplot of CD4Tcm", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(plot.title = element_text(size = 16),
                          axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
                          axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  )+ geom_smooth(method="loess",formula = y ~ x, 
                 size=1,se=FALSE, #添加置信区间             
                 color="black",              
                 linetype="solid",              
                 aes(group=1))
p33

eosinophils <- long_data2[which(long_data2$cell_type == "eosinophils"),] 
p34 = ggplot(eosinophils, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.05) +  # 设置y轴范围 
  labs(title = "Eosinophil", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))
p34

CD4Temra <- long_data2[which(long_data2$cell_type == "CD4Temra"),] 
head(CD4Temra) 
p35 = ggplot(CD4Temra, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.2) +  # 设置y轴范围 
  labs(title = "CD4Temra", x = "Group", y = "Cell Ratios") + 
  theme_minimal()+ 
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
              size=1,se=FALSE, #添加置信区间             
              color="black",              
              linetype="solid",              
              aes(group=1)) 
p35

CD8Temra <- long_data2[which(long_data2$cell_type == "CD8Temra"),] 
head(CD8Temra) 
p36 =ggplot(CD8Temra, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.02) +  # 设置y轴范围 
  labs(title = "CD8Temra", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
              size=1,se=FALSE, #添加置信区间             
              color="black",              
              linetype="solid",              
              aes(group=1))

NMonocyte <- long_data2[which(long_data2$cell_type == "NMonocyte"),] 
head(NMonocyte) 
p37 = ggplot(NMonocyte, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.0125) +  
  labs(title = "NMonocyte", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                  size=1,se=FALSE, #添加置信区间             
                  color="black",              
                  linetype="solid",              
                  aes(group=1))

M2 <- long_data2[which(long_data2$cell_type == "M2"),] 
p38 = ggplot(M2, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.015) +  # 设置y轴范围 
  labs(title = "M2", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

Tr1 <- long_data2[which(long_data2$cell_type == "Tr1"),] 
p39 =ggplot(Tr1, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.014) +  # 设置y轴范围 
  labs(title = "Tr1", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    axis.title.x = element_blank()  
  ) + geom_smooth(method="loess",formula = y ~ x, 
                            size=1,se=FALSE,      
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

regulatoryNK <- long_data2[which(long_data2$cell_type == "regulatoryNK"),] 
p40=ggplot(regulatoryNK, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.012) +  # 设置y轴范围 
  labs(title = "regulatoryNK", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

plasma <- long_data2[which(long_data2$cell_type == "plasma"),] 
p41=ggplot(plasma, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.25) +  # 设置y轴范围 
  labs(title = "PC", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

plasmablast <- long_data2[which(long_data2$cell_type == "plasmablast"),] 
p42 =ggplot(plasmablast, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.01125) +  # 设置y轴范围 
  labs(title = "PB", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

Treg <- long_data2[which(long_data2$cell_type == "Treg"),] 
p43=ggplot(Treg, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.015) +  # 设置y轴范围 
  labs(title = "Treg", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

exhaustedB <- long_data2[which(long_data2$cell_type == "exhaustedB"),] 
p44 =ggplot(exhaustedB, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.4) +  # 设置y轴范围 
  labs(title = "Bex", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

eosinophils <- long_data2[which(long_data2$cell_type == "eosinophils"),] 
p45=ggplot(eosinophils, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.05) +  # 设置y轴范围 
  labs(title = "Eosinophil", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1))

cDC1 <- long_data2[which(long_data2$cell_type == "cDC1"),] 
p46 =ggplot(cDC1, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.0113) +  # 设置y轴范围 
  labs(title = "cDC1", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

cDC2 <- long_data2[which(long_data2$cell_type == "cDC2"),] 
p47 =ggplot(cDC2, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 1, color = "black") +  # 添加散点
  ylim(0, 0.04) +  # 设置y轴范围 
  labs(title = "cDC2", x = "Group", y = "Cell Ratios") + 
  theme_minimal() + theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # 调整x轴标签字体大小和倾斜角度
    axis.title.x = element_blank()  # 如果不需要x轴标题，可以将其隐藏
  ) + geom_smooth(method="loess",formula = y ~ x, #线性拟合，二次则改为formula = y~poly(x,2)   
                            size=1,se=FALSE, #添加置信区间             
                            color="black",              
                            linetype="solid",              
                            aes(group=1)) 

combined_plot4 <- plot_grid(p28, p29, p30, p31, p32,
                            p33, p34, p35, p36, p37,
                            p38, p39, p40, p41, p42,
                            p43, p44, p45, p46, p47) 

ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/FigureSupplement.pdf",
       combined_plot4, width = 18, height = 15, device = cairo_pdf) 

### 把变化放在一起用热图表示
resultage2_normalized$group = group
resultage2_grouped_mean <- resultage2_normalized %>%
  group_by(group) %>% 
  summarize(across(everything(), mean, na.rm = TRUE)) 
dim(resultage2_grouped_mean)
resultage2_grouped_mean[1:3,]
resultage2_grouped_mean <- resultage2_grouped_mean %>% as.data.frame()
resultage2_grouped_mean <- resultage2_grouped_mean %>% column_to_rownames("group")

library(ComplexHeatmap) 
result_matrix <- as.matrix(resultage2_grouped_mean)
result_matrix_scaled <- as.data.frame(scale(resultage2_grouped_mean))
result_matrix_scaled[1:3,1:3] 
colnames(result_matrix_scaled) 
result_matrix_scaled = result_matrix_scaled[,c("CMonocyte","M0","monoDC","CD4Tnaive",
                                               "CD8Tnaive","mast_cell","ILC2","basophils","Breg",
                                               "neutrophils","CD4Tem","CD8Tem","FOB",
                                               "GC_B","gdT","ILC1",
                                               "MAIT","memoryB","MZB","NKT","Th1","Th2",
                                               "Tfh","Bnaive","Treg","NMonocyte","pDC","Tr1",
                                               "CD4Temra","CD8Temra","cDC1","cDC2","exhaustedB",
                                               "Th1/Th17","Th17","ILC3","CD8Tcm","CD4Tcm",
                                               "M2","regulatoryNK","plasmablast","eosinophils",
                                               "plasma","cytotoxicNK",
                                               "exhausted_T","IMonocyte","M1","Tc")]
result_matrix_scaled <- as.data.frame(t(result_matrix_scaled))

p5=Heatmap(as.matrix(result_matrix_scaled),
        name = "Expression",  # 热图标题
        cluster_rows=F,
        cluster_columns=F,
        show_row_names = TRUE,  # 显示行名称（细胞类型）
        show_column_names = TRUE,  # 显示列名称（group）
        col = colorRampPalette(c('navy',"white",'red'))(100),  # 颜色渐变从蓝到红
        row_title = "Cell Types",  # 行标题
        column_title = "Groups",  # 列标题
        row_dend_side = "left",    # 行的树状图显示在左侧
        column_dend_side = "top", # 列的树状图显示在上方
        heatmap_legend_param = list(title = "Expression", at = c(0, 1), labels = c("Low", "High"))
)

pdf("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/Figure3heatmap.pdf", 
    width = 4, height = 10)  # 设置 PDF 尺寸
p5                  # 绘制热图
dev.off()  

########### 
head(resultage2_normalized) 
resultage2_normalized$Bcells <- rowSums(resultage2_normalized[, c("Bnaive", "Breg", "FOB", "GC_B",
                                                                  "MZB", "memoryB", "plasma", 
                                                                  "plasmablast","exhaustedB")])
resultage2_normalized$Tcells <- rowSums(resultage2_normalized[, c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tnaive",
                                                                  "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tnaive",
                                                                  "MAIT", "NKT", "Tc", "Tfh", "Th1", "Th1/Th17",
                                                                  "Th17", "Th2", "Tr1", "Treg", "exhausted_T",
                                                                  "gdT")])  
resultage2_normalized$ILCs <- rowSums(resultage2_normalized[, c("ILC1", "ILC2", "ILC3")]) 
resultage2_normalized$NK <- rowSums(resultage2_normalized[, c("regulatoryNK","cytotoxicNK")]) 
resultage2_normalized$Macrophages <- rowSums(resultage2_normalized[, c("M0", "M1", "M2")]) 
resultage2_normalized$Granulocytes <- rowSums(resultage2_normalized[, c("neutrophils", "basophils", "eosinophils", "mast_cell")]) 
resultage2_normalized$DCs <- rowSums(resultage2_normalized[, c("pDC", "cDC1", "cDC2", "monoDC")]) 
resultage2_normalized$Monocytes <- rowSums(resultage2_normalized[, c("CMonocyte", "IMonocyte", "NMonocyte")])  

resultage_normalized3 <- resultage2_normalized[,c("Bcells","Tcells","ILCs","NK",
                                                  "Macrophages","Granulocytes","DCs","Monocytes")]
resultage_normalized3[1:3,] 
resultage_normalized3$group = resultage2_normalized$group
resultage_normalized3 <- resultage_normalized3 %>% group_by(group) %>% dplyr::summarize_each(funs(mean)) 
resultage_normalized3
View(resultage_normalized3) 
resultage_normalized3 <- resultage_normalized3 %>% as.data.frame()
resultage_long <- resultage_normalized3 %>%
  pivot_longer(cols = -group, names_to = "cell_type", values_to = "proportion")
head(resultage_long)

cell_order <- c(
  "Tcells", "Bcells", "ILCs","NK",       # 淋巴细胞优先
  "Granulocytes","Monocytes", "Macrophages", "DCs"
)
resultage_long$cell_type <- factor(
  resultage_long$cell_type,
  levels = cell_order
) 

library(RColorBrewer)
cell_colors <- brewer.pal(8, "Set2") 
resultage_long
resultage_long <- resultage_long %>%
  mutate(group = case_when(
    group == "group0to1" ~ "age0-1",
    group == "group10+" ~ "age10-20",
    group == "group20+" ~ "age20-30",
    group == "group30+50" ~ "age30-50",
    group == "group50to70" ~ "age50-70",
    group == "group70+" ~ "age70+",
    TRUE ~ group  # 保持其他值不变
  )) 

pE = ggplot(resultage_long, aes(x = group, y = proportion, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_colors) +
  labs(
    x = "Group",
    y = "Proportion",
    title = "Distribution of Immune Cell Types Across Age Groups"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box.spacing = unit(0.5, "cm")
  ) 
pE
pdf("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/age/Figure3barplot.pdf", 
    width = 6.5, height = 7)  # 设置 PDF 尺寸
p6                 # 绘制热图
dev.off() 
