getwd()
setwd("<LOCAL_R_ROOT>/Fig4/")
library(data.table) 
library(dplyr)
library(tibble)
BiocManager::install("qs")
library(qs) 
library(dplyr) 
library(tibble) 
result_TCGA <- fread(file = "<LOCAL_R_ROOT>/Fig4/TCGA_BRCA_ImmuCellAI2_marker9000/TCGA_BRCA_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
result_TCGA <- result_TCGA %>% column_to_rownames("V1") 
result_TCGA[1:3,1:3] 
dim(result_TCGA) 
which(is.na(result_TCGA)) 
result_TCGA <- result_TCGA[rowSums(is.na(result_TCGA)) == 0,] 

data= readRDS("<LOCAL_R_ROOT>/Fig4/data6690.RData") 
data[1:3,1:3] 
dim(data) 
result_TCGA <- result_TCGA[rownames(data),] 
dim(result_TCGA) 

d <- result_TCGA %>% t() %>% as.data.frame() %>% as.matrix() 
seed <- 2020 
maxK <- 6 # maximum number of clusters to try
getwd() 
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus) 

getwd() 
results <- ConsensusClusterPlus(d, distance = "euclidean",
                                maxK = maxK, reps = 1000, 
                                pItem = 0.8, pFeature = 1, 
                                title = "cluster_TCGA",
                                innerLinkage = "complete", 
                                seed = seed, plot = "pdf", 
                                clusterAlg = "pam", corUse = "complete.obs", 
                                writeTable = T, 
                                tmyPal = grDevices::colorRampPalette(c("#000004FF", 
                                                                       "#56106EFF", 
                                                                       "#BB3754FF", 
                                                                       "#F98C0AFF", 
                                                                       "#FCFFA4FF"))(64)) 
str(results) 
consensusMatrix_select <- results[[4]][["consensusMatrix"]] 
dim(consensusMatrix_select) 
consensusMatrix_select[1:5,1:5] 

library(qs) 
df <- qread("result_TCGA_select_ZN.qs") 
# write.csv(names(df), "colnames_heatmap.csv", row.names = F, col.names = F)  
df2 <- read.csv("colnames_heatmap.csv") 
names(df)[-1] <- df2$name 

expr <- df 
expr[1:3,1:3] 
df[1:3,1:3] 
cluster <- read.csv("<LOCAL_R_ROOT>/Fig4/cluster_TCGA/cluster_TCGA.k=5.consensusClass.csv") 
head(cluster) 
head(cluster44_class) 
result_TCGA[1:3,1:3] 
expr = result_TCGA
expr <- expr %>% rownames_to_column("sample") 

df <- cluster44_class %>% 
  inner_join(expr) %>% 
  column_to_rownames(names(.)[1]) 
df[1:3,1:3]
FUN <- function(gene) {
  temp <- df[, c(gene, "cluster")]
  names(temp)[1] <- "gene"
  pv <- kruskal.test(gene ~ cluster, data = temp)$p.value
  return(data.frame(gene = gene, pv = pv))
} 
l <- lapply(names(df)[-1], FUN)
res <- do.call(rbind, l)
p.value <- res$pv
res$sig.label <- ifelse(p.value < 0.001, "****",
                        ifelse(p.value < 0.005, "***",
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*", "")
                               )
                        )
)
# res$lab <- str_c(res$gene, res$sig.label)
res$lab <- res$gene


# heatmap -----------------------------------------------------------------
expr[1:3,1:3] 
head(cluster44_class) 
expr <- expr %>%
  column_to_rownames("sample") %>%
  t() %>%
  as.data.frame()
# rownames(expr) <- res$lab

# exprSet <- expr

c1_samples <- cluster44_class$sample[cluster44_class$cluster == "1"] 
c2_samples <- cluster44_class$sample[cluster44_class$cluster == "2"] 
c3_samples <- cluster44_class$sample[cluster44_class$cluster == "3"] 
c4_samples <- cluster44_class$sample[cluster44_class$cluster == "4"] 

set.seed(NULL)

c1_shuffled <- sample(c1_samples)
c2_shuffled <- sample(c2_samples)
c3_shuffled <- sample(c3_samples)
c4_shuffled <- sample(c4_samples)

desired_order <- c(c1_shuffled, c2_shuffled, c3_shuffled, c4_shuffled)

desired_order <- desired_order[desired_order %in% cluster44_class$sample]
cluster44_class <- cluster44_class[match(desired_order, cluster44_class$sample), ]

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
exprSet <- expr_df 

# expr_df <- expr[, match(cluster44_class$sample,colnames(expr))] 
# expr_df <- expr[, cluster44_class$sample]  
exprSet <- expr_df 
exprSet[1:3,1:3] 
n <- t(scale(t(exprSet))) 
n[n > 1] <- 1 
n[n < -1] <- -1 
exprSet <- n %>% as.data.frame()  

# exp <- exprSet[, cluster %>% dplyr::filter(cluster == "Cluster3") %>% .$sample]
# sel <- exp["basophils****", ] %>% as.numeric()
# sel <- which(sel > 0) 
# sel <- names(exp)[sel] 
# save(sel, file = "sel_sample.RData") 
# 
# identical(names(exprSet), cluster$sample)

exprSet[1:3,1:3] 
exprSet <- exprSet[c("BGC","M0","PB","PC","CD4Temra","Treg","pDC","CD8Trm", 
                     "Eosinophil","ILC1","Tc","Th1","Th1/17","Th17","Th2","Tr1", 
                     "Mast cell","M1","NKreg","cMo","MDSC","Neutrophil","CD4Trm","TAM",
                     "intMo","CD8Temra","Basophil", 
                     "Bex","Bnaive","Breg","CD4Tcm","CD4Tem","CD4Tn","CD8Tcm","CD8Tem", 
                     "CD8Tn","FOB","MZB","ILC2","ILC3","Langerhans","M2","MAIT", 
                     "MBC","NKT","Tex","Tfh","cDC1","cDC2","cNK","gdT","monoDC","ncMo"),] 
library(ComplexHeatmap) 
column_ha <- HeatmapAnnotation(cluster = anno_block( 
  gp = gpar(fill = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF", "#59A14FFF",
                     "#EDC948FF", "#B07AA1FF", "#FF9DA7FF", "#9C755FFF") %>% .[1:5], lwd = 2), 
  labels = c("IP1", "IP2", "IP3", "IP4"),
  labels_gp = gpar(col = "white", fontsize = 10)
)) 

Heatmap(exprSet, 
             border = TRUE,
             col = colorRampPalette(c("#3341A3", "black", "#E0DA54"))(50),
             border_gp = gpar(col = "black", lwd = 2),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             top_annotation = column_ha,
             column_split = cluster4_class$cluster,
             column_title = NULL, # <- this removes 1, 2, 3
             show_column_names = FALSE, width = ncol(exprSet) * unit(0.02, "mm"), 
             height = nrow(exprSet) * unit(3.5, "mm"), 
             heatmap_legend_param = list(title = "", title_gp = gpar(fontsize = 8), 
                                         labels_gp = gpar(fontsize = 11)), 
             row_names_gp = gpar(fontsize = 10) 
) 

pdf("heatmap_cluster_TCGA.pdf", width = 8, height = 8)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(p, heatmap_legend_side = "left", annotation_legend_side = "left", newpage = F, merge_legend = T)
popViewport()
dev.off()
