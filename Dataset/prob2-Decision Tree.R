### Single decision tree
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

## Fitting a single decision tree model
library(tree)

tree.fit <- tree(surv_ind ~.-sample_id, data)
{par(mfrow = c(1,1))
  plot(tree.fit)
  text(tree.fit,pretty = 0)
  tree.fit}


## Prunning
set.seed(3)
cv.treefit = cv.tree(tree.fit, FUN = prune.misclass)
cv.treefit

# Plot the error rate (dev) as a function of both size and k
{par(mfrow = c(1,2))
  plot(cv.treefit$size,cv.treefit$dev,type = "b")
  plot(cv.treefit$k,cv.treefit$dev,type = "b")}

##Fitting Classification Trees
# which(cv.treefit$dev == min(cv.treefit$dev)) # Index for minimum dev
prune.treefit = prune.misclass(tree.fit, best = cv.treefit$size[tail(which(cv.treefit$dev == min(cv.treefit$dev)),1)])
{plot(prune.treefit)
  text(prune.treefit,pretty = 0)}

## Cross-validation (K-folds)
# load the library
library(caret)

K <- 5 # Number of folds of the cross-validation
folds <- createFolds(1:nrow(data), K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)
error_rate <- c()
for(i in 1:K){
  fold_ind <- folds[[i]]
  data.train <- data[-fold_ind, ]
  data.test <- data[fold_ind, ]
  
  model <- tree(surv_ind ~.-sample_id, data.train)
  model.prune <- cv.tree(model, FUN = prune.misclass)
  model.best <- prune.misclass(model, best = model.prune$size[tail(which(model.prune$dev == min(model.prune$dev)),1)])
 
  tree.test <- predict(model.best, newdata=data.test, type = "class")
  survival_index <- data.test$surv_ind
  truth.table <- table(tree.test, survival_index)
  
  error_rate <- append(error_rate, 1-sum(diag(truth.table))/sum(truth.table))
}

error_rate; mean(error_rate); sd(error_rate)