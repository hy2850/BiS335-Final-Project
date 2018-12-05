### Bagging tree
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



## Fitting a bagging tree
library(randomForest)

bag.fit <- randomForest(surv_ind ~.-sample_id, data = data, mtry = 868, importance = TRUE)
{par(mfrow = c(1,1))
  plot(bag.fit)
  bag.fit}

## Cross-validation (K-folds)
# load the library
library(caret)

K <- 5 # Number of folds of the cross-validation
obs_num <- 1:nrow(data)
folds <- createFolds(obs_num, K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)

error_mean_per <- c()
tree_size <- seq(100,1000,25)

for (size in tree_size){
  error_rate <- c()
  for(i in 1:K){6
    test.ind <- folds[[i]]
    train.ind <- setdiff(obs_num, test.ind)
    
    bag.model <- randomForest(surv_ind ~.-sample_id, data = data, subset = train.ind, mtry = 868, importance = TRUE, ntree = size)
    
    bag.test <- predict(bag.model, newdata=data[test.ind,], type = "class")
    survival_index <- data[test.ind,]$surv_ind
    truth.table <- table(bag.test, survival_index)
    
    error_rate <- append(error_rate, 1-sum(diag(truth.table))/sum(truth.table))
  }
  
  error_mean_per <- append(error_mean_per, mean(error_rate))
}

tree_size; error_mean_per;

# n.tree does not seem important in reducing error after some value
# Try changing mtry - use tuneRF
# https://stackoverflow.com/questions/13956435/setting-values-for-ntree-and-mtry-for-random-forest-regression-model




## <ntree ì§ì ‘ ? •?•´ì£¼ëŠ” ê²½ìš°>
# error_rate <- c()
# for(i in 1:K){
#   test.ind <- folds[[i]]
#   train.ind <- setdiff(obs_num, test.ind)
#   
#   #  fold.ind <- folds[[i]]
#   #  data.train <- data[-fold_ind, ]
#   #  data.test <- data[fold_ind, ]
#   
#   bag.model <- randomForest(surv_ind ~.-sample_id, data = data, subset = train.ind, mtry = 868, importance = TRUE, ntree = 200)
#   
#   bag.test <- predict(bag.model, newdata=data[test.ind,], type = "class")
#   survival_index <- data[test.ind,]$surv_ind
#   truth.table <- table(bag.test, survival_index)
# 
#   error_rate <- append(error_rate, 1-sum(diag(truth.table))/sum(truth.table))
# }
# 
# error_rate; mean(error_rate); sd(error_rate)


# ## Random Forest
# error_rate <- c()
# for(i in 1:K){
#   test.ind <- folds[[i]]
#   train.ind <- setdiff(obs_num, test.ind)
# 
#   #  fold.ind <- folds[[i]]
#   #  data.train <- data[-fold_ind, ]
#   #  data.test <- data[fold_ind, ]
# 
#   rand.model <- randomForest(surv_ind ~.-sample_id, data = data, subset = train.ind, mtry = sqrt(868), importance = TRUE)
# 
#   rand.test <- predict(rand.model, newdata=data[test.ind,], type = "class")
#   survival_index <- data[test.ind,]$surv_ind
#   truth.table <- table(rand.test, survival_index)
# 
#   error_rate <- append(error_rate, 1-sum(diag(truth.table))/sum(truth.table))
# }
# 
# error_rate; mean(error_rate); sd(error_rate)