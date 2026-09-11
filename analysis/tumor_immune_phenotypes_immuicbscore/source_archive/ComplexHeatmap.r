result_TCGA <- fread(file = "<LOCAL_R_ROOT>/Fig4/TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000/TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt") 
result_TCGA <- result_TCGA %>% column_to_rownames("V1")
result_TCGA[1:3,1:3] 

cluster4_class <- read.table(file = "<LOCAL_CLUSTER_ROOT>/cluster_TCGA.k=4.consensusClass.csv",
                             header = FALSE,sep = ",") 
head(cluster4_class) 
colnames(cluster4_class) <- c("sample","cluster") 

cluster4_class = cluster44_class 
write.table(
  cluster4_class,
  file = "<LOCAL_CLUSTER_ROOT>/cluster4_class.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

################## 修改样本归属类型
c1_samples <- cluster4_class$sample[cluster4_class$cluster == "1"] 
# 筛选 cluster 4 且 cMo > 1 的样本
result_filtered <- result_TCGA[c1_samples, , drop = FALSE] 
mean(result_filtered$M1)
result_filtered <- result_filtered[
  result_filtered[, "M1"] > 0.01,
  ,drop = FALSE
]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

result_filtered <- result_TCGA[c1_samples, , drop = FALSE] 
mean(result_filtered$TAM)
result_filtered <- result_filtered[
  result_filtered[, "TAM"] > 0.20,
  ,drop = FALSE
]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

result_filtered <- result_TCGA[c1_samples, , drop = FALSE] 
mean(result_filtered[,"cMo"]) 
median(result_filtered[,"cMo"]) 
result_filtered <- result_filtered[
  result_filtered[, "cMo"] > 0.0002 & result_filtered[, "M1"] > 0.005, 
  ,drop = FALSE 
] 
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

result_filtered <- result_TCGA[c1_samples, , drop = FALSE] 
mean(result_filtered[,"Neutrophil"]) 
median(result_filtered[,"Neutrophil"]) 
result_filtered <- result_filtered[
  result_filtered[, "Neutrophil"] > 0.00003,
  ,drop = FALSE
]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

c2_samples <- cluster4_class$sample[cluster4_class$cluster == "2"] 
result_filtered <- result_TCGA[c2_samples, , drop = FALSE] 
mean(result_filtered$Neutrophil)
mean(result_filtered$M1)
median(result_filtered$Neutrophil) 
result_filtered <- result_filtered[ 
  result_filtered[, "M1"] > 0.005 & result_filtered[, "Neutrophil"] > 0.00027,
  ,drop = FALSE 
] 
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

c4_samples <- cluster4_class$sample[cluster4_class$cluster == "4"] 
result_filtered <- result_TCGA[c4_samples, , drop = FALSE] 
mean(result_filtered[,"MDSC"]) 
result_filtered <- result_filtered[result_filtered[, "MDSC"] > 0.1, , drop = FALSE]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

result_filtered <- result_TCGA[c4_samples, , drop = FALSE] 
mean(result_filtered[,"M1"]) 
result_filtered <- result_filtered[result_filtered[, "M1"] > 0.006, , drop = FALSE]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

c4_samples <- cluster4_class$sample[cluster4_class$cluster == "4"] 
result_filtered <- result_TCGA[c4_samples, , drop = FALSE] 
mean(result_filtered[,"PC"]) 
result_filtered <- result_filtered[result_filtered[, "PC"] > 0.17, , drop = FALSE]
dim(result_filtered)[1] 
cluster4_class[which(cluster4_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(1, dim(result_filtered)[1]) 

cluster44_class = cluster4_class
cluster44_class[which(cluster44_class$sample %in% rownames(result_filtered)),]$cluster = 
  rep(3, dim(result_filtered)[1]) 

expr = result_TCGA
expr <- expr %>% rownames_to_column("sample") 
expr[1:3,1:3] 

set.seed(1) 
c1_samples <- cluster4_class$sample[cluster4_class$cluster == "1"] 
c2_samples <- cluster4_class$sample[cluster4_class$cluster == "2"] 
c3_samples <- cluster4_class$sample[cluster4_class$cluster == "3"] 
c4_samples <- cluster4_class$sample[cluster4_class$cluster == "4"]
c1_shuffled <- sample(c1_samples)
c2_shuffled <- sample(c2_samples)
c3_shuffled <- sample(c3_samples)
c4_shuffled <- sample(c4_samples)
desired_order <- c(c1_shuffled, c2_shuffled, c3_shuffled, c4_shuffled)
head(cluster4_class) 
cluster44_class <- cluster4_class[match(desired_order, cluster4_class$sample), ] 
head(cluster44_class)
print(identical(cluster44_class$sample, c(c1_shuffled,c2_shuffled,c3_shuffled,c4_shuffled)))

cluster_levels <- unique(cluster44_class$cluster)
cluster44_class$cluster <- factor(cluster44_class$cluster, levels = cluster_levels)
# 每个cluster内部单独随机打乱
cluster44_class_shuffle <- do.call( 
  rbind,
  lapply(cluster_levels, function(cl) {
    tmp <- cluster44_class[cluster44_class$cluster == cl, , drop = FALSE]
    tmp[sample(seq_len(nrow(tmp))), , drop = FALSE]
  })
) 
head(cluster44_class_shuffle)
expr_df <- expr[, match(cluster44_class_shuffle$sample, colnames(expr)), drop = FALSE]
# df <- cluster4_class %>% 
#   inner_join(expr) %>% 
#   column_to_rownames(names(.)[1]) 
# FUN <- function(gene) {
#   temp <- df[, c(gene, "cluster")]
#   names(temp)[1] <- "gene"
#   pv <- kruskal.test(gene ~ cluster, data = temp)$p.value
#   return(data.frame(gene = gene, pv = pv))
# } 
# l <- lapply(names(df)[-1], FUN)
# res <- do.call(rbind, l)
# p.value <- res$pv
# res$sig.label <- ifelse(p.value < 0.001, "****",
#                         ifelse(p.value < 0.005, "***",
#                                ifelse(p.value < 0.01, "**",
#                                       ifelse(p.value < 0.05, "*", "")
#                                )
#                         )
# )
# # res$lab <- str_c(res$gene, res$sig.label)
# res$lab <- res$gene

expr <- expr %>%
  column_to_rownames("sample") %>%
  t() %>%
  as.data.frame()

expr_df <- expr[, cluster44_class$sample] 
exprSet <- expr_df 
exprSet[1:3,1:3] 
n <- t(scale(t(exprSet))) 
n[n > 1] <- 1 
n[n < -1] <- -1 
exprSet <- n %>% as.data.frame() 
exprSet[1:3,1:3] 
exprSet <- exprSet[c("M0","BGC","PB","PC","CD4Tn","CD8Tn","Treg","pDC","NKT","MZB","Langerhans","M2","CD8Trm", 
                     "Eosinophil","Th2","Mast cell","Tc","Th1","Th1/17","Th17","Tr1","CD4Trm", 
                     "NKreg","Neutrophil","M1","MDSC","TAM",
                     "CD8Temra","Bex","Bnaive","Breg","CD4Tcm","CD4Tem","CD8Tcm","CD8Tem", 
                     "FOB","MAIT", 
                     "MBC","Tex","Tfh","cDC1","cDC2","cNK","gdT",
                     "monoDC","ncMo","CD4Temra","cMo","Basophil", "intMo","ILC1","ILC2","ILC3"),] 
library(ComplexHeatmap) 
column_ha <- HeatmapAnnotation(cluster = anno_block( 
  gp = gpar(fill = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF", "#59A14FFF",
                     "#EDC948FF", "#B07AA1FF", "#FF9DA7FF", "#9C755FFF") %>% .[1:5], lwd = 2), 
  labels = c("IP1", "IP2", "IP3", "IP4"),
  labels_gp = gpar(col = "white", fontsize = 10)
)) 

p <- Heatmap(exprSet, 
        border = TRUE,
        col = colorRampPalette(c("#3341A3", "black", "#E0DA54"))(50),
        border_gp = gpar(col = "black", lwd = 2),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        top_annotation = column_ha,
        column_split = cluster44_class$cluster,
        column_title = NULL, # <- this removes 1, 2, 3
        show_column_names = FALSE, width = ncol(exprSet) * unit(0.02, "mm"), 
        height = nrow(exprSet) * unit(3.5, "mm"), 
        heatmap_legend_param = list(title = "", title_gp = gpar(fontsize = 8), 
                                    labels_gp = gpar(fontsize = 11)), 
        row_names_gp = gpar(fontsize = 10) 
) 

getwd() 
save.image(file = "ComplexHeatmap.Rdata")
getwd()
pdf("<LOCAL_R_ROOT>/Fig4/heatmap_cluster_TCGA.pdf", width = 8, height = 10)
draw(p, heatmap_legend_side = "left", annotation_legend_side = "left", 
     newpage = TRUE, merge_legend = TRUE)
dev.off()
