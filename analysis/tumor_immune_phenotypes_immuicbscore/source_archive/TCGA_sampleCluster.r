getwd()
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/case_TCGA/") 
load("ConsensusMatrix.Rdata")
load("annoConsensusMatrix.Rdata") 
load("exprSet_matrix.Rdata") 
result_TCGA <- read.xlsx("<LOCAL_LEGACY_PROJECT_ROOT>/backup0106/learning/myresult/ResultDeconvolutionTCGAtumor.xlsx",
                         sheet = 1) 
result_TCGA <- fread(file = "<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/case_TCGA/TCGA_ImmuCellAI2_state_fraction.txt")    
result_TCGA[1:3,1:5] 
result_TCGA <- result_TCGA %>% column_to_rownames("V1")
dim(result_TCGA)

library(qs)
library(data.table)
library(dplyr) 

data <- qread("result_TCGA_cancer.qs") 
data[1:3,1:3] 
dim(data)
result_TCGA <- result_TCGA[rownames(data),]

data <- data %>% t() %>% as.data.frame()

boxplot(data[,1:20],main="before") 
boxplot(exprSet[,1:20],main="after") 
class(exprSet)
exprSet <- exprSet %>% t() %>% as.data.frame()

exprSet[1:3,1:3]
exprSet <- as.matrix(exprSet)
exprSet <- apply(exprSet, 2, as.numeric)

dim(exprSet)
exprSet <- exprSet[complete.cases(exprSet), ]  

exprSet <- scale(exprSet)
View(exprSet) 
class(exprSet)
dim(exprSet)
rownames(exprSet) <- rownames(data)

library(ConsensusClusterPlus) 
class(exprSet)
exprSet[1:3,1:3]
exprSet <- as.data.frame(t(exprSet))
dist_matrix <- dist(exprSet, method = "canberra")
View(dist_matrix)
print(dist_matrix)

matrix(exprSet)
result_TCGA_select[1:3,1:3]
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/case_TCGA/") 
clinical <- read.delim("ann_TCGA.csv",sep = ",") 
clinical[1:3,] 
dim(clinical) 
unique(clinical$tissue) 
clinical <- clinical %>% dplyr::select(-X) 

result_TCGA_select <- result_TCGA_select %>% t() %>% as.data.frame() 
clinical <- clinical[which(clinical$sample %in% colnames(result_TCGA_select)),] 
rownames(clinical) <- clinical$sample 
cancer_clinical <- clinical[which(clinical$type2 == "tumor"),] 
dim(cancer_clinical) 
cancer_clinical[1:3,] 
cancer_clinical <- cancer_clinical[,-1] 
unique(cancer_clinical$tissue) 
clinical_COAD <- clinical[which(clinical$project_id == "TCGA-COAD"),]
clinical_LIHC <- clinical[which(clinical$project_id == "TCGA-LIHC"),]
clinical_SKCM <- cancer_clinical[which(cancer_clinical$tissue == "SKCM"),]
dim(clinical_SKCM)
clinical_COAD[1:3,1:3] 
carcinoma_clinicalTumor  <- cancer_clinical[which(cancer_clinical$tissue %in% c("SKCM")),]
dim(carcinoma_clinicalTumor)
carcinoma_clinicalTumor[1:3,1:3]
dim(result_TCGA) 
result_TCGA[1:3,1:3]
result_TCGA_select <- result_TCGA[,rownames(carcinoma_clinical)] 
dim(result_TCGA_select) 
result_TCGA_select[1:3,1:3]
# result_TCGA_select <- result_TCGA_select[-53,] 
rownames(result_TCGA_select)
result_TCGA_select <- as.matrix(result_TCGA_select) 
dim(result_TCGA_select) 
result_TCGA_select <- result_TCGA_select %>% t() %>% as.data.frame() 
result_TCGA_select[1:3,1:3] 
result_liver[,1:3]
rownames(result_TCGA_select) 
# result_TCGA_select <- result_TCGA_select[c("M1","M2","regulatoryNK","cytotoxicNK","exhausted_T","CD8Trm","CD8Tnaive",
#                                            "CD4Tnaive","CD4Trm","exhaustedB","Bnaive","CD8Tem","plasma","plasmablast",
#                                            "cDC1","cDC2","Treg","TAM"),] 
getwd()
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/case_TCGA/cluster_TCGA_Solid/")
setwd("<LOCAL_LEGACY_PROJECT_ROOT>/CaseStudy/case_TCGA/")
qsave(result_TCGA_select, file = "result_TCGA_select.qs")
result_TCGA_select[1:3,1:3]
dim(result_TCGA_select)
dim(result_TCGA)
cluster3Cells2 = result_TCGA[,cluster3Cells2$V1]
cluster3Cells2[1:3,1:3]
dim(cluster3Cells2)
dim(cluster3Cells)
getwd()
result_TCGA <- result_TCGA %>% t() %>% as.data.frame()
result_TCGA[1:3,1:3]
result_TCGA <- result_TCGA %>% drop_na()

getwd()
results <- ConsensusClusterPlus(as.matrix(result_TCGA), 
                                maxK=5, 
                                reps=1000, 
                                pItem=0.8, 
                                pFeature=1, 
                                title="cluster4", 
                                clusterAlg="pam", 
                                distance="euclidean", 
                                seed=123456, 
                                writeTable = T,
                                plot="pdf",innerLinkage = "complete", 
                                tmyPal = grDevices::colorRampPalette(c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF"))(64)) 
dim(results)

getwd() 
load("Ensembl2Symbol.RData") 
consensusMatrix_select <- results[[2]][["consensusMatrix"]] 
consensusMatrix_select[1:5,1:5] 
dim(consensusMatrix_select) 
colnames(consensusMatrix_select) <- colnames(cluster3Cells2) 
rownames(consensusMatrix_select) <- colnames(cluster3Cells2) 

ConsensusMatrix <- data.frame(results[[2]][["consensusMatrix"]]) 
ConsensusMatrix <- ConsensusMatrix[results[[2]]$consensusTree$order,
                                   results[[2]]$consensusTree$order]
#创建注释列数据框
annCol <- data.frame(results=paste0("Cluster", 
                                    results[[2]][['consensusClass']] 
                                    [results[[2]]$consensusTree$order]), 
                     row.names=colnames(ConsensusMatrix)) 
annColors<-list(results=c("Cluster1"="#db6968", 
                          "Cluster2"="#4d97cd",
                          "Cluster3"="#99cbeb",
                          "Cluster4"="#459943", 
                          "Cluster4"="#56106EFF")) 
BiocManager::install("pheatmap") 
library(pheatmap) 
ConsensusMatrix[1:5,1:5] 
dim(ConsensusMatrix)
xx <- as.data.frame(matrix(nrow = 1919,ncol = 0)) 
xx$Celltype <- rownames(consensusMatrix_select)
head(xx)
xx$number <- 1:1919
rownames(xx) <- 1:1919
xx <- xx[rownames(ConsensusMatrix),]
head(ConsensusMatrix)
dim(ConsensusMatrix)
rownames(ConsensusMatrix) <- xx$Celltype
colnames(ConsensusMatrix) <- xx$Celltype
ConsensusMatrix[1:3,1:3]
dim(ConsensusMatrix) 
Heatmap <- pheatmap(ConsensusMatrix,
                    color=colorRampPalette((c("white","steelblue")))(100),
                    cluster_cols=FALSE,
                    cluster_rows=FALSE,
                    clustering_distance_cols="correlation",
                    clustering_method="average",
                    border_color=NA,
                    annotation_col=annCol,
                    annotation_colors=annColors,
                    show_colnames=TRUE,
                    show_rownames=TRUE) 

##########  TCGA数据分期
setwd("<LOCAL_LEGACY_AUX_ROOT>/Tcell_states/1.data/TCGA/GDCdata/clinical/") 
clinical <- read.delim("clinical.tsv") 
head(clinical);dim(clinical) 
clinical$case_submitter_id  
result_TCGA_select[1:3,1:3] 
result_TCGA_select <- result_TCGA_select %>% t() %>% as.data.frame() 
result_TCGA_select$sampleID <- rownames(result_TCGA_select) %>% str_extract(.,".*(?=-)") 
head(result_TCGA_select$sampleID) 
dim(result_TCGA_select)
length(intersect(result_TCGA_select$sampleID,clinical$case_submitter_id)) 
clinical <- clinical[which(clinical$case_submitter_id %in% result_TCGA_select$sampleID),]
dim(clinical) 
head(clinical$case_submitter_id)                                                                                                                                                                                                                      
clinical[1:3,1:3] 
clinical <- clinical[,c("case_id","case_submitter_id","project_id","age_at_index","ajcc_pathologic_stage",
                        "ajcc_pathologic_t","ajcc_pathologic_m","ajcc_pathologic_n")]
clinical_unique <- distinct(clinical) 
clinical_unique[1:3,];dim(clinical_unique) 

