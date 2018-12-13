library(MASS)


clin <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/clinical.rds")
gex <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/expression.rds")
mut <- readRDS("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/mutation.rds")
genes <- read.csv("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/gex_anova_result.csv")
mutat <- read.csv("C:/Users/YuKyeong/Desktop/Lecture/BiS335/final_PJ/Finalterm-Project/mut_chisq_result.csv")
surv_ind <- clin$survival_index
genes <- genes[order(genes[,2]),]
mutat <- mutat[order(mutat[,2]),]


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


# PCA 
data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
rownames(data) <- common_id
data <- cbind(data,t(gex))
data <- cbind(data, binary)

k <- prcomp(data[,-1], center = T, scale = T)
data <- data.frame(surv_ind = surv_ind[match(common_id, clin$sample_id)])
rownames(data) <- common_id
data <- cbind(data, k$x)
data <- rbind(data, data[data$surv_ind==2, ]) # repeat surv_ind = 2 data

# SVM
# forward n backward hybrid
library(e1071)
library(MASS)
library(permute)

candi <- c(2:481)
selected <- c(1)
candi <- setdiff(candi, selected)
temp2 <- temp <- c(0)

while(1==1){
  rate <- c()
  for (j in candi){
    acc.test <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,append(selected, j)]
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
      test <- data[indices,]
      svmfit <- svm(surv_ind~., data = new_data, kernel = "radial", subset = -indices)
      pred.test <- predict(svmfit, newdata = test)
      acc.test[i] <- mean(test$surv_ind == pred.test)
      # tune <- tune(svm, surv_ind~., data = new_data[-indices, ], kernel = "linear", ranges = list(cost = c(0.01, 0.1, 1, 10), gamma = c(0.5, 1, 2, 4))) 
      # pred.test <- predict(tune$best.model, newdata = test)
      # acc.test[i] <- mean(test$surv_ind == pred.test)
    }
    rate <- append(rate, mean(acc.test))
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
    acc.test <- c()
    random <- shuffle(nrow(data))
    new_data <- data[,selected[-j]]
    for (i in 1:5){
      indices <- random[(((i-1)*round(length(random)/5))+1):(i*round(length(random)/5))]
      test <- data[indices,]
      svmfit <- svm(surv_ind~., data = new_data, kernel = "radial", subset = -indices)
      pred.test <- predict(svmfit, newdata = test)
      acc.test[i] <- mean(test$surv_ind == pred.test)
      # tune <- tune(svm, surv_ind~., data = new_data[-indices, ], kernel = "linear", ranges = list(cost = c(0.01, 0.1, 1, 10), gamma = c(0.5, 1, 2, 4))) 
      # pred.test <- predict(tune$best.model, newdata = test)
      # acc.test[i] <- mean(test$surv_ind == pred.test)
    }
    rate2 <- append(rate2, mean(acc.test))
  }
  if (temp[length(temp)] < max(rate2)){
    print(max(rate2))    
    candi <- sort(append(candi, selected[-1][which.max(rate)]))
    selected <- selected[-(which.max(rate2)+1)]
    temp <- append(temp, max(rate2))
  }else if (mark == 1){break}
  print(selected)
}

