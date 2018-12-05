### Tree Boosting
rm(list = ls()) # Clear all workspace
setwd("C:/BiS335 Final Project/Dataset")

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
# Boosting (Using built-in k-fold CV method and gbm.perf)
library(gbm)

# set.seed(1)
# train <- sample(1:nrow(data), nrow(data)*0.75)
# data.train <- data[train,]
# data.test <- data[-train,]
# 
# # "Bernouli" distribution -> For 0,1 response / Use "multinomial" distribution for 1~4 response classification?
# # predict.gbm -> returns probability? to be classified into each class defined in gbm model with "multinomial" dist
# # Q. How to control the parameters? n.trees, interaction.dept, and shrinkage?
# #   -> 1. n.trees can be found out using "gbm.perf"
# boost.fit <- gbm(surv_ind ~.-sample_id, data = data.train, distribution = "multinomial", n.trees = 500, shrinkage = 0.001, interaction.dept = 4, cv.folds = 5)
# tree.num <- gbm.perf(boost.fit, method = "cv")
# 
# boost.pred <- predict(boost.fit, newdata = data.test, n.trees = tree.num, type = "response")
# boost.test <-apply(boost.pred, 1, which.max) # select class with highest probability
# 
# survival_index <- data.test$surv_ind
# truth.table <- table(boost.test, survival_index)
# 
# error_rate <- 1-sum(diag(truth.table))/sum(truth.table)
# error_rate


#########################################################################################################
## Cross-validation (K-folds)
# load the library
library(caret)

K <- 5 # Number of folds of the cross-validation
obs_num <- 1:nrow(data)
folds <- createFolds(obs_num, K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)

tree_size <- 500
error_rate <- c()

for(i in 1:K){
   fold_ind <- folds[[i]]
   data.train <- data[-fold_ind, ]
   data.test <- data[fold_ind, ]

   boost.model <- gbm(surv_ind ~.-sample_id, data = data.train, distribution = "multinomial", n.trees = tree_size, shrinkage = 0.001, interaction.dept = 4)

   tree.num <- gbm.perf(boost.fit)
   boost.pred <- predict(boost.fit, newdata = data.test, n.trees = tree.num, type = "response")

#   boost.pred <- predict(boost.model, newdata=data.test, n.trees = tree_size)
   boost.test <-apply(boost.pred, 1, which.max) # select class with highest probability

   survival_index <- data.test$surv_ind
   truth.table <- table(boost.test, survival_index)

   error_rate <- append(error_rate, 1-sum(diag(truth.table))/sum(truth.table))
}

error_rate; mean(error_rate); sd(error_rate)

# Current best : Boosting size 500, shrinkage 0.001, interaction 4, gbm.perf OOB