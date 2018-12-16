library(MASS)
# install.packages("permute")
library(permute)
library(e1071)
library(gbm)
library(caret)

clin <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/clinical.rds") # data load
mut <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/mutation.rds")
gex <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/expression.rds")
genes <- c("ZFHX3", "ZFP36", "FCGR1A", "DHRS9", "FGF18", "RPL36A", "SLC4A2", "TEX264", "SERPINF2", "OAZ2",
           "PSAP", "ENAM", "IL26", "GPR4", "PPP1R10")
genes2 <- c("MAP3K1", "IL27RA", "ZFP36", "UQCRB", "LOC388272", "STAU2", "SLC2A13", "FLJ40298", "PSMD14",
            "KRT34", "LOC197322", "SYNJ2", "WBSCR17", "POLH", "M6PR", "MGMT")
mutat2 <- c("FABP6", "NPAP1", "MAPKBP1")

surv_ind <- clin$survival_index

gex <- na.omit(gex)
common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex)) # reduce data
gex1 <- gex[which(rownames(gex) %in%  genes), match( common_id,colnames(gex))]
gex2 <- gex[which(rownames(gex) %in%  genes2), match( common_id,colnames(gex))]
mut <- mut[which(mut[,1] %in% common_id),c(1,2)]
mut1 <- mut[which(mut[,2] %in% mutat2), ]
mut2 <- mut[which(mut[,2] %in% mutat2), ]

# make binary matrix for mutation
binary <- matrix(data = 0, nrow = length(common_id), ncol = length(mutat2))
colnames(binary) <- paste0(mutat2, 'mut')
for (i in 1:length(mutat2)){
  binary[match(mut[,1][mut[,2]==mutat2[i]], common_id), i] <- 1
}

# data frame
data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
rownames(data) <- common_id
data1 <- cbind(data,t(gex1))
data2 <- cbind(data,t(gex2))
data2 <- cbind(data2, binary)

# divide test and train data
test_indices <- shuffle(nrow(data))[1:round(nrow(data)*0.2)]
test1 <- data1[test_indices, ]
data1 <- data1[-test_indices, ]
test2 <- data2[test_indices, ]
data2 <- data2[-test_indices, ]
error.test <- c() # save error rate for each model

#----------------lda----------------
lda.fit <- lda(surv_ind ~., data = data1)
pred.test <- predict(lda.fit, test1, type = "response")
acc.test <- mean(test1$surv_ind == pred.test$class)
print(paste0("LDA error rate : ", 1-acc.test))

lda.fit <- lda(surv_ind ~., data = data2[,-c(18:20)])
pred.test <- predict(lda.fit, test2, type = "response")
acc.test <- mean(test2$surv_ind == pred.test$class)
print(paste0("LDA error rate : ", 1-acc.test))

#---------------qda-----------------

# qda.fit <- qda(surv_ind ~., data = data1)
# pred.test <- predict(qda.fit, test1, type = "response")
# acc.test <- mean(test1$surv_ind == pred.test$class)
# print(paste0("LDA error rate : ", 1-acc.test))
# 
# qda.fit <- qda(surv_ind ~., data = data2[,-c(18:20)])
# pred.test <- predict(qda.fit, test2, type = "response")
# acc.test <- mean(test2$surv_ind == pred.test$class)
# print(paste0("LDA error rate : ", 1-acc.test))

#--------------svm------------------

svm.fit <- svm(surv_ind ~., data = data1, kernel = "polynomial")
pred.test <- predict(svm.fit, newdata = test1)
acc.test <- mean(test1$surv_ind == pred.test)
print(paste0("SVM error rate : ", 1-acc.test))

svm.fit <- svm(surv_ind ~., data = data2, kernel = "polynomial")
pred.test <- predict(svm.fit, newdata = test2)
acc.test <- mean(test1$surv_ind == pred.test)
print(paste0("SVM error rate : ", 1-acc.test))

#------------ boosting -------------
tree_size <- 2000 # Default : 100
shrinkage.control = 0.001 # Default : 0.1
interaction.dept.control = c(6) # Default : 1

acc.test <- c()

boost.fit <- gbm(surv_ind ~., data = data1, distribution = "multinomial", 
                 n.trees = tree_size, shrinkage = shrinkage.control, 
                 interaction.dept = interaction.dept.control, cv.folds = 5, 
                 bag.fraction = 1, n.minobsinnode = 18)
tree.num <- gbm.perf(boost.fit, method = "cv")

boost.pred <- predict(boost.fit, test1, n.trees = tree.num, type = "response")
boost.test <- apply(boost.pred, 1, which.max)

acc.test <- mean(boost.test == test1$surv_ind)
print(paste0("boosting error rate : ", 1-acc.test))

boost.fit <- gbm(surv_ind ~., data = data2, distribution = "multinomial", 
                 n.trees = tree_size, shrinkage = shrinkage.control, 
                 interaction.dept = interaction.dept.control, cv.folds = 5, 
                 bag.fraction = 1, n.minobsinnode = 18)
tree.num <- gbm.perf(boost.fit, method = "cv")

boost.pred <- predict(boost.fit, test2, n.trees = tree.num, type = "response")
boost.test <- apply(boost.pred, 1, which.max)

acc.test <- mean(boost.test == test2$surv_ind)
print(paste0("boosting error rate : ", 1-acc.test))