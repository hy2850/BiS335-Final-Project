### Single decision tree
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
set.seed(1)
train <- sample(1:nrow(data), nrow(data)*0.7)
data.train <- data[train,]
data.test <- data[-train,]


## Cross-validation (K-folds)
# load the library
library(caret)
library(tree)

error_rate <- c()

K <- 5 # Number of folds of the cross-validation
folds <- createFolds(1:nrow(data.train), K) # Requires 'caret' package; Creates K-folds for testing (Using rest of the data for training)
for(i in 1:K){
  fold_ind <- folds[[i]]
  data.training <- data.train[-fold_ind, ]
  data.validation <- data.train[fold_ind, ]
  
  
  model <- tree(surv_ind ~.-sample_id, data.training)
  model.prune <- cv.tree(model, FUN = prune.misclass)
  model.best <- prune.misclass(model, best = model.prune$size[tail(which(model.prune$dev == min(model.prune$dev)),1)])
  
  tree.test_cv <- predict(model.best, newdata=data.validation, type = "class")
  survival_index_cv <- data.validation$surv_ind
  truth.table_cv <- table(tree.test_cv, survival_index_cv)
  
  error_rate <- append(error_rate, 1-mean(survival_index_cv==tree.test_cv))
}

error_rate; mean(error_rate); sd(error_rate)

#########################################################################################################
## Building real model and testing with data.test
## Fitting a single decision tree model

tree.fit <- tree(surv_ind ~.-sample_id, data.train)
{par(mfrow = c(1,1))
  plot(tree.fit)
  text(tree.fit,pretty = 0)
  #tree.fit
  }


## Prunning
set.seed(3)
cv.treefit = cv.tree(tree.fit, FUN = prune.misclass)
#cv.treefit

# Plot the error rate (dev) as a function of both size and k
{par(mfrow = c(1,2))
  plot(cv.treefit$size,cv.treefit$dev,type = "b")
  plot(cv.treefit$k,cv.treefit$dev,type = "b")}

##Fitting Classification Trees
# which(cv.treefit$dev == min(cv.treefit$dev)) # Index for minimum dev
prune.treefit = prune.misclass(tree.fit, best = cv.treefit$size[tail(which(cv.treefit$dev == min(cv.treefit$dev)),1)])
{plot(prune.treefit)
  text(prune.treefit,pretty = 0)}

tree.test <- predict(prune.treefit, newdata=data.test, type = "class")
survival_index <- data.test$surv_ind

error_rate_final <- 1-mean(survival_index==tree.test); error_rate_final

##########################################################################
## Feature selection process
library(MASS)
library(permute)

candi <- c(3:870) # predicator candidate

selected <- c(1,2) # selected predictor 
temp <- c(0) # current max accuracy 

while(1==1){
  ## 1. Forward stepwise selection
  rate <- c() # accuracy for each 5-fold CV
  for (j in candi){
    acc.test <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,append(selected, j)] # add j th predicator data
    
    # 5-Fold CV    
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
      test <- data[indices,]
      
      ####################### Change this to apply to other models ##########################
      model <- tree(surv_ind ~. -sample_id, data = new_data, subset = -indices)
      #model.prune <- cv.tree(model, FUN = prune.misclass)
      #model.best <- prune.misclass(model, best = model.prune$size[tail(which(model.prune$dev == min(model.prune$dev)),1)])
      
      pred.test <- predict(model, test, type = "class")
      #######################################################################################
      
      # conf.test <- table(test$surv_ind, pred.test$class);
      acc.test[i] <- mean(test$surv_ind == pred.test)
    }
    
    rate <- append(rate, mean(acc.test))
  }
  
  if (temp[length(temp)] < max(rate)){ # If newly selected feature increases the max accuracy of the model, select that feature
    print(max(rate))
    selected <- append(selected, candi[which.max(rate)])
    selected <- sort(selected)
    temp <- append(temp, max(rate))
    candi <- candi[-which.max(rate)]
    mark <- 0 # This 'mark' boolean data is needed for stop condition of outer-most while-loop
  }else{mark <- 1} # mark = 1 indicates that no more features are selected to increase the model accuracy
  
  if (length(selected) <= 2) {next} # skip if selected has 2 predicators ()
  
  ## 2. Backward stepwise selection (almost same as forward selection, but some parameters are different)
  rate2 <- c() # accuracy for each 5-fold CV (in backward selection)
  for (j in 2:length(selected[-1])){ # Excluded recently added feature, because removing it right after it is being added is pointless
    acc.test <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,selected[-j]] # Exclude one feature from selected
    
    # 5-Fold CV   
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
      test <- data[indices,]
      
      ####################### Change this to apply to other models ##########################
      model <- tree(surv_ind ~. -sample_id, data = new_data, subset = -indices)
      #model.prune <- cv.tree(model, FUN = prune.misclass)
      #model.best <- prune.misclass(model, best = model.prune$size[tail(which(model.prune$dev == min(model.prune$dev)),1)])
      
      pred.test <- predict(model, test, type = "class")
      #######################################################################################
      
      # conf.test <- table(test$surv_ind, pred.test$class);
      acc.test[i] <- mean(test$surv_ind == pred.test)
    }
    
    rate2 <- append(rate2, mean(acc.test))
  }
  if (temp[length(temp)] < max(rate2)){ # If deleting certain feature increases the max accuracy of the model, delete that feature completely
    print(max(rate2))    
    candi <- sort(append(candi, selected[-1][which.max(rate)]))
    selected <- selected[-(which.max(rate2)+1)]
    temp <- append(temp, max(rate2))
  }else if (mark == 1){break} # If no feature increases model accuracy when deleted and also the mark is 1, finish the program
  print(selected)
}