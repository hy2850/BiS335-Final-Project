### Tree Boosting
rm(list = ls()) # Clear all workspace
setwd("C:/Users/hy285/Desktop/BiS335-Final-Project/Dataset")

## Environment settings
# setwd to directory containing below rds files!
clin <- readRDS("clinical.rds")
mut <- readRDS("mutation.rds")
gex <- readRDS("expression.rds")
genes <- read.csv("gex_anova_result.csv")
mutat <- read.csv("mut_chisq_result.csv")
surv_ind <- clin$survival_index
genes <- genes[order(genes[,2]),]
mutat <- mutat[order(mutat[,2]),]

library(MASS)

common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex)) # reduce data
gex <- gex[which(rownames(gex) %in%  genes[,1]), match( common_id,colnames(gex))]
gex <- na.omit(gex)
mut <- mut[which(mut[,1] %in% common_id),c(1,2)]
mut <- mut[which(mut[,2] %in% mutat[,1]), ]

# how do we handle na in gex?

# make binary matrix for mutation
binary <- matrix(data = 0, nrow = length(common_id), ncol = length(mutat[,1]))
colnames(binary) <- paste0(mutat[,1], 'mut')
for (i in 1:length(mutat[,1])){
  binary[match(mut[,1][mut[,2]==mutat[i,1]], common_id), i] <- 1
}

# data frame
data <- data.frame(sample_id = common_id)
data$surv_ind <- surv_ind[match(common_id, clin$sample_id)]
data <- cbind(data,t(gex))
data <- cbind(data, binary)

# Fixing bad characters in the column names
column_names <- matrix(data = colnames(data), ncol = length(colnames(data)))
colnames(data)<- as.vector(apply(column_names, 2, FUN = make.names))

#########################################################################################################
## Data separation (Train and Test)
library(gbm)

set.seed(4321)
train <- sample(1:nrow(data), nrow(data)*0.9)
data.train <- data[train,]
data.test <- data[-train,]

#########################################################################################################
## Cross-validation (K-folds)
# load the library
library(caret)
library(gbm)

error_rate <- c()

# Control gbm parameters here
tree_size <- 500 # Default : 100
shrinkage.control = 0.001 # Default : 0.1
interaction.dept.control = 6 # Default : 1

# size 500, shrinkage 0.001, interaction 2 or 4 : 0.4903

# K <- 5 # Number of folds of the cross-validation
# folds <- createFolds(1:nrow(data.train), K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)
# for(i in 1:K){
#    fold_ind <- folds[[i]]
#    data.training <- data.train[-fold_ind, ]
#    data.validation <- data.train[fold_ind, ]
#    
#    boost.model_cv <- gbm(surv_ind ~.-sample_id, data = data.training, distribution = "multinomial", n.trees = tree_size, shrinkage = shrinkage.control, interaction.depth = interaction.dept.control, bag.fraction = 0.5, n.minobsinnode = 10)
#    
#    tree.num <- gbm.perf(boost.model_cv)
#    boost.pred_cv <- predict(boost.model_cv, newdata = data.validation, n.trees = tree.num, type = "response")
#    
#    #boost.pred <- predict(boost.model, newdata=data.test, n.trees = tree_size)
#    boost.test_cv <-apply(boost.pred_cv, 1, which.max) # select class with highest probability
#    
#    survival_index_cv <- data.validation$surv_ind
#    
#    truth.table_cv <- table(boost.test_cv, survival_index_cv) # For debugging
#    error_rate <- append(error_rate, 1-mean(survival_index_cv==boost.test_cv))
# }
# 
# error_rate; mean(error_rate); sd(error_rate)

#########################################################################################################
## Building real model and testing with data.test
## Boosting (Using built-in k-fold CV method and gbm.perf)

# "Bernouli" distribution -> For 0,1 response / Use "multinomial" distribution for 1~4 response classification?
# predict.gbm -> returns probability? to be classified into each class defined in gbm model with "multinomial" dist
# Q. How to control the parameters? n.trees, interaction.depht, and shrinkage?
#   -> 1. n.trees can be found out using "gbm.perf"
boost.fit <- gbm(surv_ind ~.-sample_id, data = data.train, distribution = "multinomial", n.trees = tree_size, shrinkage = shrinkage.control, interaction.dept = interaction.dept.control, cv.folds = 5, bag.fraction = 0.8, n.minobsinnode = 20)
tree.num <- gbm.perf(boost.fit, method = "cv")

boost.pred <- predict(boost.fit, newdata = data.test, n.trees = tree.num, type = "response")
boost.test <-apply(boost.pred, 1, which.max) # select class with highest probability

survival_index <- data.test$surv_ind
truth.table <- table(boost.test, survival_index)

error_rate_final <- 1-mean(survival_index==boost.test); error_rate_final
importance <- summary.gbm(boost.fit, plotit=TRUE)


#########################################################################################################
# Lowest error rate recorded : 0.3958333
# Parameters used :
# train <- sample(1:nrow(data), nrow(data)*0.9)
# 
# tree_size <- 500
# shrinkage.control = 0.001
# interaction.dept.control = 6
# 
# bag.fraction = 0.8, n.minobsinnode = 20