# install.packages("factoextra") + NbClust
library(factoextra) 
library(NbClust)

clin <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/clinical.rds") # data load
mut <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/mutation.rds")
gex <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/expression.rds")
# genes <- read.csv("C:/Users/Lee/Desktop/Finalterm-Project/gex_anova_result.csv")
# mutat <- read.csv("C:/Users/Lee/Desktop/Finalterm-Project/mut_chisq_result.csv")
surv_ind <- clin$survival_index
# genes <- genes[order(genes[,2]),]
# mutat <- mutat[order(mutat[,2]),]

gex <- na.omit(gex)
common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex)) # reduce data
gex <- gex[, match( common_id,colnames(gex))]
# mut <- mut[which(mut[,1] %in% common_id),c(1,2)]
# mut <- mut[which(mut[,2] %in% mutat[,1]), ]

# data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
# data <- cbind(data,t())

set.seed(2) # kmean clustering
pca_gex <- prcomp(t(gex), scale = T) # pca 
pve <- 100*pca_gex$sdev^2/sum(pca_gex$sdev^2)

km.out <- kmeans(pca_gex$x, 2, nstart = 20)

par(c(1,2)) # visualization silhouetter and etc. measure for optimal number of clusters
fviz_nbclust(pca_gex$x, kmeans, method = "silhouette") 
fviz_nbclust(pca_gex$x, kmeans, method = "wss")
fviz_nbclust(pca_gex$x, kmeans, method = "gap_stat")

# NbSluster function
NbClust(data = pca_gex$x, distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans")
