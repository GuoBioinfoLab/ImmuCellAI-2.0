getwd()
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/excel_2/") 

result_Tan <- read.xlsx("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/Result_20250418.xlsx",
                        sheet = 1) %>% column_to_rownames("X1") 
result_Tan[1:3,1:3] 
dim(result_Tan) 

library(readxl) 
library(dplyr)
library(purrr)
library(tibble)
folder_path <- "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/excel_2"
file_list <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
file_list 

process_file <- function(file_path) {
  tryCatch({
    df <- read_excel(file_path) %>% as.data.frame()
    required_cols <- c("Run", "BioProject", "Disease")
    present_cols <- intersect(required_cols, colnames(df))
    if (length(present_cols) < 3) {
      missing_cols <- setdiff(required_cols, present_cols)
      warning(paste("文件", basename(file_path), "缺少列:", paste(missing_cols, collapse = ", ")))
      return(data.frame(Run = character(), BioProject = character(), Disease = character()))
    }
    result <- as.data.frame(df) %>% dplyr::select(all_of(required_cols))
    result$source_file <- basename(file_path)
    return(result)
  }, error = function(e) {
    message(paste("处理文件", basename(file_path), "时出错:", e$message))
    return(data.frame(Run = character(), BioProject = character(), Disease = character()))
  }) 
} 
combined_data <- map_dfr(file_list, process_file)
cat("共处理", length(file_list), "个文件\n")
cat("合并后的数据有", nrow(combined_data), "行，", ncol(combined_data), "列\n")
combined_data[1:3,]
unique(combined_data$Disease)

# 肺结核
Tuberculosis_data = combined_data[which(combined_data$Disease %in% c("Tuberculosis")),] 
dim(Tuberculosis_data) 

PRJNA31975 <- read_excel("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/excel_2/PRJNA638653.xlsx") %>% 
  as.data.frame() %>% dplyr::select(c("Run","BioProject","Disease","tb_status"))  
PRJNA395234 <- read_excel("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/excel_2/PRJNA395234.xlsx") %>% 
  as.data.frame() %>% dplyr::select(c("Run","BioProject","Disease","condition"))  
PRJNA31975[1:3,] 
PRJNA395234[1:3,] 
PRJNA638653[1:3,] 
unique(PRJNA31975$tb_status)
colnames(PRJNA395234)[4] = "tb_status" 
PRJNA_combined <- rbind(PRJNA31975, PRJNA395234)  
result_combined <- result_Tan[, PRJNA_combined$Run]

table(PRJNA_combined$tb_status)
dim(result_combined) 
result_combined[1:3,1:3] 
result_combined = result_combined %>% t() %>% as.data.frame() 
dim(result_combined) 
result_combined$status = PRJNA_combined$tb_status 
dim(result_combined)
unique(result_combined$status) 
result_combined[which(result_combined$status == "TB"),]$status = "activeTB"
result_combined$status <- gsub("latent TB infection", "LTBI", result_combined$status)
unique(PRJNA_combined$tb_status)
PRJNA_combined[which(PRJNA_combined$tb_status %in% c("LTBI","latent TB infection")),]$tb_status <- "Latent TB"
PRJNA_combined[which(PRJNA_combined$tb_status %in% c("activeTB","TB")),]$tb_status <- "Active TB"
PRJNA_combined[1:3,]
dim(PRJNA_combined)
getwd()
write_csv(PRJNA_combined, file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/PRJNA_combined.csv")


result_combined[1:3,1:3]
result_combined$exchange = result_combined$plasma
result_combined$plasma = result_combined$plasmablast
result_combined$plasmablast = result_combined$plasma
dim(result_combined)
result_combined <- result_combined[,-53] 
data_long <- pivot_longer(result_combined, cols = -status, names_to = "CellType", values_to = "Abundance")
head(data_long)  
unique(data_long$CellType)

selected_cells <- c("cDC1", "Tc", "plasma", "ILC1")  
data_long <- data_long[data_long$CellType %in% selected_cells, ]

getwd() 
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/") 
save(data_long,file = "data_long.Rdata") 
load("data_long.Rdata") 
library(ggsci) 
library(ggpubr) 
data_long$CellType <- factor( 
  data_long$CellType, 
  levels = c("cDC1", "Tc", "plasma", "ILC1") 
) 
p0 <- ggplot(data_long, aes(x = CellType, y = Abundance)) +
  stat_boxplot(aes(colour = status),
               geom = 'errorbar',
               position = position_dodge(width = 0.8),
               width = 0.4) +
  geom_boxplot(aes(colour = status),
               width = 0.6,
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(data = data_long, 
    aes(fill = status),
              shape = 21,
              size = 3,
              position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~ CellType, nrow = 1, scales = 'free') +
  stat_compare_means(
    aes(group = status),
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = "top",
    size = 8,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )) + theme_bw() +
  scale_fill_npg() +
  scale_color_manual(values = c('black', 'black')) +
  labs(x = NULL) +
  theme(axis.text = element_text(color = 'black', size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(size = 14, color = 'black'),  # ← 分面标题设置
        strip.background = element_blank(),  # 分面标题背景框
        axis.title = element_text(size = 14, color = 'black', face = 'bold'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) 
p0 
getwd() 
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/") 
ggsave("Tuberculosis.pdf", p0, width = 18, height = 6, 
       device = cairo_pdf, family = "Arial", limitsize = FALSE) 

############### AIDS ################ 
AIDs_data = combined_data[which(combined_data$Disease %in% c("AIDS")),] 
dim(AIDs_data) 
unique(AIDs_data$BioProject) 
PRJNA683803 <- read_excel("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/excel_2/PRJNA683803.xlsx") %>% 
  as.data.frame() %>% dplyr::select(c("Run","BioProject",
                                      "Disease","outcome","treatment")) 
head(PRJNA683803) 
unique(PRJNA683803$outcome) 
write_csv(PRJNA683803, file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/PRJNA683803_combined.csv")

result_combined <- result_Tan[, PRJNA683803$Run] 
dim(result_combined) 
result_combined[1:3,1:3]  
result_combined = result_combined %>% t() %>% as.data.frame() 
dim(result_combined) 
result_combined$status = PRJNA683803$outcome  
result_combined[1:3,] 
colnames(result_combined)
data_long <- pivot_longer(result_combined, cols = -status, names_to = "CellType", values_to = "Abundance") 
head(data_long) 
data_long$status <- factor(
  data_long$status, 
  levels = c("died", "survived", "died-IRIS", "survived-IRIS") 
) 

data_long_filtered <- data_long %>%
  group_by(CellType, status) %>%
  mutate(Q1 = quantile(Abundance, 0.25),
         Q3 = quantile(Abundance, 0.75),
         IQR = Q3 - Q1,
         lower = Q1 - 1.5 * IQR,
         upper = Q3 + 1.5 * IQR) %>%
  filter(Abundance >= lower & Abundance <= upper)
library(ggpubr)
library(dplyr)

my_comparisons <- list(
c("died", "survived"),
c("died-IRIS", "survived-IRIS")
)
stat.test <- compare_means(
Abundance ~ status,
group.by = "CellType",
data = data_long_filtered,
method = "wilcox.test",
comparisons = my_comparisons
)
# 添加 y 位置
y_positions <- data_long_filtered %>%
group_by(CellType) %>%
summarise(max_y = max(Abundance, na.rm = TRUE)) %>%
mutate(y.position = max_y + 0.05 * max_y) %>%
dplyr::select(CellType, y.position)
# 合并并补充坐标信息
stat.test <- stat.test %>%
left_join(y_positions, by = "CellType") %>%
mutate(
xmin = group1,
xmax = group2
)
stat.test
# 绘图
stat.test = stat.test %>%
  filter((group1 == "died" & group2 == "survived") |
           (group1 == "died-IRIS" & group2 == "survived-IRIS"))
head(data_long_filtered) 
unique(data_long_filtered$CellType) 
# data_long_main <- data_long_filtered 
data_long_main <- data_long_filtered[which(data_long_filtered$CellType %in% 
                                             c("cytotoxicNK","Th1", "exhausted_T", 
                                               "Tc","CD8Temra",
                                               "M2")),] 
# stat.test_main <- stat.test
stat.test_main <- stat.test[which(stat.test$CellType %in% c("cytotoxicNK","Th1", 
                                                            "exhausted_T","Tc",
                                                            "CD8Temra","M2")),]

data_long_main$CellType <- factor(data_long_main$CellType, 
                                  levels = c("cytotoxicNK", "Th1", "exhausted_T", 
                                             "Tc", "CD8Temra", "M2")) 
p = ggplot(data_long_main, aes(x = status, y = Abundance)) +
stat_boxplot(aes(colour = status),
geom = 'errorbar',
position = position_dodge(width = 0.8),
width = 0.4) +
geom_boxplot(aes(colour = status),
width = 0.6,
outlier.shape = NA,
position = position_dodge(width = 0.8)) +
geom_jitter(aes(fill = status, colour = status),
shape = 21,
position = position_jitterdodge(jitter.width = 0.4),
alpha = 0.7) +
facet_wrap(~ CellType, nrow = 4, scales = 'free') +
stat_pvalue_manual(
stat.test_main,
label = "p.signif",
tip.length = 0, 
hide.ns = FALSE
) + 
theme_bw() + 
scale_fill_npg() + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 14, face = "bold"),  # 增大分面标题字体
        axis.title = element_text(size = 14),  # 增大坐标轴标题字体（可选）
        legend.text = element_text(size = 12) 
  ) 
p
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/AIDs_mainFigure.pdf",
       p, width = 50, height = 14, device = cairo_pdf, limitsize = FALSE)

## 附图 
data_long_supplement <- data_long_filtered[which(!data_long_filtered$CellType %in% 
                                             c("basophils","neutrophils",
                                               "CMonocyte","cytotoxicNK",
                                               "Tc","Th1","Th17","Tr1")),] 
stat.test_supplement <- stat.test[which(!stat.test$CellType %in% c("basophils","neutrophils",
                                                            "CMonocyte","cytotoxicNK",
                                                            "Tc","Th1","Th17","Tr1")),]
p2 = ggplot(data_long_supplement, aes(x = status, y = Abundance)) +
  stat_boxplot(aes(colour = status),
               geom = 'errorbar',
               position = position_dodge(width = 0.8),
               width = 0.4) +
  geom_boxplot(aes(colour = status),
               width = 0.6,
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = status),
              shape = 21,
              position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~ CellType, scales = 'free') +
  stat_pvalue_manual(
    stat.test_supplement,
    label = "p.signif",
    tip.length = 0, 
    hide.ns = FALSE
  ) + 
  theme_bw() + 
  scale_fill_npg() + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 14, face = "bold"),  # 增大分面标题字体
        axis.title = element_text(size = 14),  # 增大坐标轴标题字体（可选）
        legend.text = element_text(size = 12) 
  ) 
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/AIDs_supplementFigure.pdf",
       p2,width = 20, height = 20, device = cairo_pdf)

# early ART VS deferred ART
result_combined[1:3,1:3] 
result_combined = result_combined %>% t() %>% as.data.frame() 
dim(result_combined) 
head(PRJNA683803)  
unique(PRJNA683803$outcome)
PRJNA683803[which(PRJNA683803$outcome %in% c("died-IRIS")),]$outcome <- "Died-IRIS"
PRJNA683803[which(PRJNA683803$outcome %in% c("survived-IRIS")),]$outcome <- "Survived-IRIS"
PRJNA683803[which(PRJNA683803$outcome %in% c("survived")),]$outcome <- "Survived"
PRJNA683803[which(PRJNA683803$outcome %in% c("died")),]$outcome <- "Died"
PRJNA683803[which(PRJNA683803$treatment %in% c("Early")),]$treatment <- "Early ART"
PRJNA683803[which(PRJNA683803$treatment %in% c("deferred")),]$treatment <- "Deferred ART"
write_csv(PRJNA683803, file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/PRJNA683803_combined.csv")


result_combined$status = PRJNA683803$treatment  
result_combined[1:3,] 
unique(PRJNA683803$treatment) 
data_long2 <- pivot_longer(result_combined, cols = -status, names_to = "CellType", values_to = "Abundance") 
head(data_long2) 
data_long2$status <- factor(
  data_long2$status, 
  levels = c("Early", "deferred")  # 强制指定顺序
) 
data_long2$Abundance <- as.numeric(data_long2$Abundance) 
data_long_filtered2 <- data_long2 %>%
  group_by(CellType, status) %>%
  mutate(Q1 = quantile(Abundance, 0.25),
         Q3 = quantile(Abundance, 0.75),
         IQR = Q3 - Q1,
         lower = Q1 - 1.5 * IQR,
         upper = Q3 + 1.5 * IQR) %>% 
  filter(Abundance >= lower & Abundance <= upper) 

library(ggpubr) 
library(dplyr) 
my_comparisons <- list(
  c("Early", "deferred")
)
stat.test <- compare_means(
  Abundance ~ status,
  group.by = "CellType",
  data = data_long_filtered2,
  method = "wilcox.test",
  comparisons = my_comparisons
)
y_positions <- data_long_filtered2 %>%
  group_by(CellType) %>%
  summarise(max_y = max(Abundance, na.rm = TRUE)) %>%
  mutate(y.position = max_y + 0.05 * max_y) %>%
  dplyr::select(CellType, y.position)
stat.test <- stat.test %>%
  left_join(y_positions, by = "CellType") %>%
  mutate(
    xmin = group1,
    xmax = group2
  )
stat.test
# 绘图
data_long_supplement2 <- data_long_filtered2[which(data_long_filtered2$CellType %in% 
                                             c("cytotoxicNK","gdT","exhausted_T",
                                               "Th1")),] 
stat.test_main2 <- stat.test[which(stat.test$CellType %in% c("cytotoxicNK","gdT","exhausted_T",
                                                             "Th1")),] 
save(data_long_supplement2,data_long_supplement, data_long_main, file = "AIDs_Analysis.Rdata")
p3 = ggplot(data_long_supplement2, aes(x = status, y = Abundance)) + 
  stat_boxplot(aes(colour = status),
               geom = 'errorbar',
               position = position_dodge(width = 0.8),
               width = 0.4) +
  geom_boxplot(aes(colour = status),
               width = 0.6,
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = status, colour = status),
              shape = 21,
              position = position_jitterdodge(jitter.width = 0.4),
              alpha = 0.7) +
  facet_wrap(~ CellType, ncol=4,scales = 'free') +
  stat_pvalue_manual(
    stat.test_main2,
    label = "p.signif",
    tip.length = 0, 
    hide.ns = FALSE
  ) + 
  theme_bw() + 
  scale_fill_npg() + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 14, face = "bold"),  # 增大分面标题字体
        axis.title = element_text(size = 14),  # 增大坐标轴标题字体（可选）
        legend.text = element_text(size = 12) 
  ) 
p3
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/AIDs_ARTFigure.pdf",
       p3,width = 20, height = 50, device = cairo_pdf, limitsize = FALSE) 

# 时间点折线图
head(PRJNA683803) 
head(result_combined) 
result_combined[1:3,1:3]
dim(result_combined) 
result_combined <- result_combined[,-52] 
result_combined$grouping = PRJNA683803$grouping 

result_combined %>% 
  mutate(group = case_when( 
    (grouping %in% c("ED-I_0","DS-I_0","ES-I_0","ES-N_0","DD-N_0","ED-N_0","DS-N_0")) ~ "D0" ,
    (grouping %in% c("ED-I_1","DS-I_1","ES-I_1","ES-N_1","DD-N_1","ED-N_1","DS-N_1"))~ "D1",
    (grouping %in% c("ED-I_4i","DS-I_4i","ES-I_4i","ES-N_4","DS-N_4"))~ "D4",
    (grouping %in% c("DS-I_8","ES-I_8","ES-N_8","DS-N_8"))~ "D8"
  )) -> result_combined 
result_combined$outcome = PRJNA683803$outcome 
head(result_combined)  
result_combined <- result_combined %>% dplyr::select(-grouping) 
result_combined <- result_combined %>%
  mutate(group_outcome = paste(group, outcome, sep = "_")) 
result_combined2 = result_combined %>% dplyr::select(-group,-outcome)

result <- result_combined2 %>% 
  group_by(group_outcome) %>% dplyr::summarize_each(funs(mean))
result <- result %>% as.data.frame()
result[1:3,1:3]
dim(result)
result_long <- pivot_longer(result, cols = -group_outcome, names_to = "CellType", values_to = "Abundance")
head(result_long)
unique(result_long$group_outcome)
result_long$group = result_long$group_outcome %>% str_extract(.,".*(?=_)") 
result_long$status = result_long$group_outcome %>% str_extract(.,"(?<=_).*") 
head(result_long) 
result_long <- result_long[which(result_long$status %in% c("died","survived")),]

result_long <- result_long[which(result_long$CellType %in% 
                                c("Tc","Th1","Treg","exhausted_T")),] 
result_long$CellType <- factor(result_long$CellType, 
                                  levels = c("Tc","Th1","Treg","exhausted_T")) 

phy.cols <- c("#A3C8DC","#EA5D2D") 
p4 = ggplot(result_long, aes(x = group, y = Abundance, group = status, color = status)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~ CellType, ncol= 4, scales = "free_y") + 
  labs(x = "Time", y = "Ratios", title = "Cell Type Ratios by Status") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),  # 设置面板背景色为白色
        plot.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)  
  ) + scale_colour_manual(values =phy.cols) 
p4 
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/AIDs_zhexianFigure.pdf",
       p4, width = 10, height = 6, device = cairo_pdf) 

p0
p
p3
p4

getwd() 
save(data_long,data_long_main,data_long_supplement,data_long_supplement2,
     stat.test_main,stat.test_supplement,stat.test_main2,result_long,
     p0,p,p2,p3,p4, file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/Figure6.Rdata")

library(cowplot) 
left_column <- plot_grid(p0, p3, p4, 
                         ncol = 1,        # 垂直排列
                         align = "v",     # 垂直对齐
                         axis = "l",      # 左侧对齐
                         rel_heights = c(1, 1, 1)) 
final_plot <- plot_grid(left_column, p, 
                        ncol = 2,        # 左右并排
                        align = "h",     # 水平对齐
                        rel_widths = c(3, 2))  # 调整左右宽度比例

# 显示或保存图形
final_plot 
ggsave("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/Tan/Figure6.pdf",
       final_plot, width = 19, height = 20, device = cairo_pdf) 


