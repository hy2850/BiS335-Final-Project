library(MASS)
# install.packages("permute")
library(permute)
library(e1071)
library(gbm)
library(caret)

clin <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/clinical.rds") # data load
mut <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/mutation.rds")
gex <- readRDS("C:/Users/Lee/Desktop/Finalterm-Project/expression.rds")
genes <- read.csv("C:/Users/Lee/Desktop/Finalterm-Project/gex_anova_result.csv")
mutat <- read.csv("C:/Users/Lee/Desktop/Finalterm-Project/mut_chisq_result.csv")
surv_ind <- clin$survival_index
genes <- genes[order(genes[,2]),]
mutat <- mutat[order(mutat[,2]),]

gex <- na.omit(gex)
common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex)) # reduce data
gex <- gex[which(rownames(gex) %in%  genes[,1]), match( common_id,colnames(gex))]
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
data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
rownames(data) <- common_id
data <- cbind(data,t(gex))
data <- cbind(data, binary)

# divide test and train data
test_indices <- shuffle(nrow(data))[1:round(nrow(data)*0.25)]
test <- data[test_indices, ]
data <- data[-test_indices, ]
error.test <- c() # save error rate for each model

#--------------------------lda----------------------------
candi <- c(2:869) # predicator candidate
selected <- c(1) # 1 = surv_ind label, selected predictor
candi <- setdiff(candi, selected) # Remove selected feature from candi
temp <- c(0) # current max accuracy 

while(1==1){
  ## 1. Forward stepwise selection
  rate <- c() # accuracy for each 5-fold CV
  for (j in candi){
    acc.val <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,append(selected, j)] # add j th predicator data

    # 5-Fold CV    
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):min(length(random), (i*round(length(random)/5)))]
      val <- data[indices,]

      ####################### Change this to apply to other models ##########################
      lda.fit <- lda(surv_ind ~., data = new_data, subset = -indices)
      pred.val <- predict(lda.fit, val, type = "response")
      #######################################################################################
      
      # conf.val <- table(val$surv_ind, pred.val$class);
      acc.val[i] <- mean(val$surv_ind == pred.val$class)
    }
    
    rate <- append(rate, mean(acc.val))
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
    acc.val <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,selected[-j]] # Exclude one feature from selected
    
    # 5-Fold CV   
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):min(length(random), (i*round(length(random)/5)))]
      val <- data[indices,]

      ####################### Change this to apply to other models ##########################
      lda.fit <- lda(surv_ind ~., data = new_data, subset = -indices)
      pred.val <- predict(lda.fit, val, type = "response")
      #######################################################################################
      
      # conf.val <- table(val$surv_ind, pred.val$class);
      acc.val[i] <- mean(val$surv_ind == pred.val$class)
    }
    
    rate2 <- append(rate2, mean(acc.val))
  }
  if (temp[length(temp)] < max(rate2)){ # If deleting certain feature increases the max accuracy of the model, delete that feature completely
    print(max(rate2))    
    candi <- sort(append(candi, selected[-1][which.max(rate)]))
    selected <- selected[-(which.max(rate2)+1)]
    temp <- append(temp, max(rate2))
  }else if (mark == 1){break} # If no feature increases model accuracy when deleted and also the mark is 1, finish the program
  print(selected[-1])
}

# estimate test error
new_data <- data[,selected]
lda.fit <- lda(surv_ind ~., data = data)
pred.test <- predict(lda.fit, test, type = "response")
acc.test <- mean(test$surv_ind == pred.test$class)
print(paste0("LDA error rate : ", 1-acc.test))
error.test[1] <- 1-acc.test

#------------------qda---------------------

# data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
# rownames(data) <- common_id
# data <- cbind(data,t(gex))
# data <- cbind(data, binary)
# 
# k <- prcomp(data[,-1], center = T, scale = T)
# data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
# rownames(data) <- common_id
# data <- cbind(data,k$x)
# data <- rbind(data, data[data$surv_ind==2, ])
# 
# candi <- c(2:480)
# selected <- c(1,3)
# candi <- setdiff(candi, selected)
# temp2 <- temp <- c(0)
# 
# while(1==1){
#   rate <- c()
#   for (j in candi){
#     acc.test <- c()
#     random <- shuffle(nrow(data))
#     new_data <- data[,append(selected, j)]
#     for (i in 1:5){
#       indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
#       test <- data[indices,]
#       qda.fit <- qda(surv_ind ~., data = new_data, subset = -indices)
#       pred.test <- predict(qda.fit, test, type = "response")
#       # conf.test <- table(test$surv_ind, pred.test$class);
#       acc.test[i] <- mean(test$surv_ind == pred.test$class)
#     }
#     rate <- append(rate, mean(acc.test))
#   }
#   if (temp[length(temp)] < max(rate)){
#     print(max(rate))
#     selected <- append(selected, candi[which.max(rate)])
#     selected <- sort(selected)
#     temp <- append(temp, max(rate))
#     candi <- candi[-which.max(rate)]
#     mark <- 0
#   }else{mark <- 1}
#   
#   if (length(selected) <= 2) {next}
#   rate2 <- c()
#   for (j in 2:length(selected[-1])){
#     acc.test <- c()
#     random <- shuffle(nrow(data))
#     new_data <- data[,selected[-j]]
#     for (i in 1:5){
#       indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
#       test <- data[indices,]
#       qda.fit <- qda(surv_ind ~., data = new_data, subset = -indices)
#       pred.test <- predict(qda.fit, test, type = "response")
#       # conf.test <- table(test$surv_ind, pred.test$class);
#       acc.test[i] <- mean(test$surv_ind == pred.test$class)
#     }
#     rate2 <- append(rate2, mean(acc.test))
#   }
#   if (temp[length(temp)] < max(rate2)){
#     print(max(rate2))    
#     candi <- sort(append(candi, selected[-1][which.max(rate)]))
#     selected <- selected[-(which.max(rate2)+1)]
#     temp <- append(temp, max(rate2))
#   }else if (mark == 1){break}
#   print(selected)
# }
# # lda.fit2 <- qda(surv_ind ~ .-sample_id, data = data)
# # d2 <- table(predict(lda.fit2)$'class', data$surv_ind); d2

#---------------------------------------------svm--------------------------------------------
candi <- c(2:869)
selected <- c(1)
candi <- setdiff(candi, selected)
temp <- c(0)

while(1==1){
  rate <- c()
  for (j in candi){
    acc.val <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,append(selected, j)]
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):min(length(random), (i*round(length(random)/5)))]
      val <- data[indices,]
      svmfit <- svm(surv_ind~., data = new_data, kernel = "polynomial", subset = -indices)
      pred.val <- predict(svmfit, newdata = val)
      acc.val[i] <- mean(val$surv_ind == pred.val)
      # tune <- tune(svm, surv_ind~., data = new_data[-indices, ], kernel = "linear", ranges = list(cost = c(0.01, 0.1, 1, 10), gamma = c(0.5, 1, 2, 4))) 
      # pred.val <- predict(tune$best.model, newdata = val)
      # acc.val[i] <- mean(val$surv_ind == pred.val)
    }
    rate <- append(rate, mean(acc.val))
  }
  if (temp[length(temp)] < max(rate)){
    print(max(rate))
    selected <- append(selected, candi[which.max(rate)])
    selected <- sort(selected)
    temp <- append(temp, max(rate))
    candi <- candi[-which.max(rate)]
    mark <- 0
  }else{mark <- 1}
  
  if (length(selected) <= 2) {next}
  rate2 <- c()
  for (j in 2:length(selected[-1])){
    acc.val <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,selected[-j]]
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):min(length(random), (i*round(length(random)/5)))]
      val <- data[indices,]
      svmfit <- svm(surv_ind~., data = new_data, kernel = "polynomial", subset = -indices)
      pred.val <- predict(svmfit, newdata = val)
      acc.val[i] <- mean(val$surv_ind == pred.val)
      # tune <- tune(svm, surv_ind~., data = new_data[-indices, ], kernel = "linear", ranges = list(cost = c(0.01, 0.1, 1, 10), gamma = c(0.5, 1, 2, 4))) 
      # pred.val <- predict(tune$best.model, newdata = val)
      # acc.val[i] <- mean(val$surv_ind == pred.val)
    }
    rate2 <- append(rate2, mean(acc.val))
  }
  if (temp[length(temp)] < max(rate2)){
    print(max(rate2))    
    candi <- sort(append(candi, selected[-1][which.max(rate)]))
    selected <- selected[-(which.max(rate2)+1)]
    temp <- append(temp, max(rate2))
  }else if (mark == 1){break}
  print(selected[-1])
}

new_data <- data[,selected]
svmfit <- svm(surv_ind ~., data = data, kernel = "linear")
pred.test <- predict(svmfit, test, type = "response")
acc.test <- mean(test$surv_ind == pred.test)
print(paste0("SVM error rate : ", 1-acc.test))
error.test[2] <- 1-acc.test

#------------------boosting---------------------

# data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
# rownames(data) <- common_id
# data <- cbind(data,t(gex))
# data <- cbind(data, binary)

# k <- prcomp(data[,-1], center = T, scale = T)
# data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
# rownames(data) <- common_id
# data <- cbind(data,k$x)
# data <- rbind(data, data[data$surv_ind==2, ])

# selected <- c(1:869)

# Control gbm parameters here
tree_size <- 2800 # Default : 100
shrinkage.control = 0.001 # Default : 0.1
interaction.dept.control = c(6) # Default : 1

acc.test <- c()
optimal_trees<- c()

boost.fit <- gbm(surv_ind ~., data = data, distribution = "multinomial", 
                 n.trees = tree_size, shrinkage = shrinkage.control, 
                 interaction.dept = interaction.dept.control, cv.folds = 5, 
                 bag.fraction = 0.8, n.minobsinnode = 20)
tree.num <- gbm.perf(boost.fit, method = "cv")
# optimal_trees[i] <- which.min(boost.model_cv$valid.error)

boost.pred <- predict(boost.fit, test, n.trees = tree.num, type = "response")
boost.test <- apply(boost.pred, 1, which.max)

acc.test <- mean(boost.test == test$surv_ind)
print(paste0("boosting error rate : ", 1-acc.test))
error.test[3] <- 1 - acc.test
importance <- summary.gbm(boost.fit, plotit=TRUE)

# -----------------------------plot----------------------------

barplot(error.test, legend.text = c("LDA", "SVM", "Boosting"), ylab = "Error rate", col = c(1,2,3), ylim = c(0,1))
