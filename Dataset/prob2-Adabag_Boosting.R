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

set.seed(1)
train <- sample(1:nrow(data), nrow(data)*0.7)
data.train <- data[train,]
data.test <- data[-train,]

#########################################################################################################
## Cross-validation (K-folds)
# load the library
library(caret)

error_rate <- c()

# Control gbm parameters here
tree_size <- 2500 # Default : 100
shrinkage.control = 0.01 # Default : 0.1
interaction.dept.control = 4 # Default : 1

K <- 5 # Number of folds of the cross-validation
folds <- createFolds(1:nrow(data.train), K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)
for(i in 1:K){
  fold_ind <- folds[[i]]
  data.training <- data.train[-fold_ind, ]
  data.validation <- data.train[fold_ind, ]
  
  
  adabag.boost <- boosting(surv_ind ~.-sample_id, data = data.training, boos = FALSE, control=rpart.control(maxdepth=3))
  
  adabag.predict <- predict.boosting(adabag.boost, newdata = data.validation)
  #adabag.predict$confusion
  #adabag.predict$error
  
  error_rate <- append(error_rate, adabag.predict$error)
}

error_rate; mean(error_rate); sd(error_rate)

#########################################################################################################
## Building real model and testing with data.test
## Boosting in 'adabag' package
# adabag.boost <- boosting(surv_ind ~.-sample_id, data = data.train, boos = FALSE, control=rpart.control(maxdepth=3))
# #adabag.boost <- boosting.cv(surv_ind ~.-sample_id, data = data.train, v = 5, control = rpart.control(cp=0.01), par = TRUE)
# 
# adabag.predict <- predict.boosting(adabag.boost, newdata = data.test)
# adabag.predict$confusion
# adabag.predict$error # This returns same result as two below codes
# 
# adabag.test <-apply(adabag.predict$prob, 1, which.max) # select class with highest probability
# adabag.error_rate <- 1-mean(data.test$surv_ind==adabag.test); adabag.error_rate