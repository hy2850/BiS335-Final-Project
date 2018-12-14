clin <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/clinical.rds")
gex <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/expression.rds")
mut <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/mutation.rds")
genes <- read.csv("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/gex_anova_result.csv")
mutat <- read.csv("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/mut_chisq_result.csv")
surv_ind <- clin$survival_index
genes <- genes[order(genes[,2]),]
mutat <- mutat[order(mutat[,2]),]

gex <- na.omit(gex)
common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex)) # reduce data
gex <- gex[, match( common_id,colnames(gex))]

set.seed(2)
pca_gex <- prcomp(t(gex), scale = T) # pca 

km.out <- kmeans(pca_gex$x, 4, nstart = 20)

# visualization silhouetter and etc. measure for optimal number of clusters
fviz_nbclust(pca_gex$x, kmeans, method = "silhouette") 
fviz_nbclust(pca_gex$x, kmeans, method = "wss")
# ---------- Hierarchial Clustering ----------

library(factoextra) 
library(NbClust)

fviz_nbclust(pca_gex$x, hcut, method = "silhouette") 
# avg. sil. width = 0.1 at k = 3

library("fpc")

hc.complete = hclust(dist(data[, -1]), method="complete")
plot(hc.complete, main="Complete Linkage", xlab="", sub="", cex=.9)
hc.cluster = cutree(hc.complete, 3)

table(hc.cluster,as.numeric(data$surv_ind))
# 112 12 129 55
# 40 3 36 55
# 22 0 3 13

library(clues)
adjustedRand(hc.cluster,as.numeric(data$surv_ind))
# Rand : 0.52058977


# ---------- K means ----------

km.out <- kmeans(pca_gex$x, 3, nstart = 20)

table(km.out$cluster, as.numeric(data$surv_ind))
# 31 2 31 28
# 110 12 123 68
# 33 1 14 27

adjustedRand(km.out$cluster, as.numeric(data$surv_ind))
# Rand : 0.51083855