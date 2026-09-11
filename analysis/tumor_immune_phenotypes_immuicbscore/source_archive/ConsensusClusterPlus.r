(load("TLS.RData"))
genes <- TLS
expr <- expr[genes, ]

library(ConsensusClusterPlus)
library(tibble)
library(qs)
library(dplyr)
d <- expr %>% as.matrix()
seed <- 2020
maxK <- 10 # maximum number of clusters to try
results <- ConsensusClusterPlus(d, distance = "euclidean",
                                maxK = maxK, reps = 1000, pItem = 0.8, pFeature = 1, title = "cluster_TCGA",
                                innerLinkage = "complete", seed = seed, plot = "pdf", 
                                clusterAlg = "pam", corUse = "complete.obs", writeTable = T, 
                                tmyPal = grDevices::colorRampPalette(c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF"))(64), 
                                legendPal = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627","#023E8A","#9D4EDD", "#223D6C","#D20A13")
) # Note that we implement consensus clustering with innerLinkage="complete".
# save(results, file = "ConsensusCluster_ssGSEA_TCGA.RData") 

results[[6]][["consensusMatrix"]][1:5,1:5]
consensusMatrix_select <- results[[6]][["consensusMatrix"]] 
consensusMatrix_select[1:3,1:3]
colnames(consensusMatrix_select)<-colnames(exprSet)
rownames(consensusMatrix_select)<-colnames(exprSet)

# We advise against using innerLinkage="average" which is the default value in this package as average linkage is not robust to outliers.
############## PAC implementation ##############
Kvec <- 3:maxK
x1 <- 0.1
x2 <- 0.9  # threshold defining the intermediate sub-interval
PAC <- rep(NA, length(Kvec))  
names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
for (i in Kvec) {
  M <- results[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
} # end for i# The optimal K
optK <- Kvec[which.min(PAC)]
optK

