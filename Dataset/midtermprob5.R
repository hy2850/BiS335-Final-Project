# put the appropriate path name in your computer
clin <- readRDS("C:/BiS335-Final-Project/Dataset/midterm/clin.Rds")
gex <- readRDS("C:/BiS335-Final-Project/Dataset/midterm/expression.rds")
mut <- readRDS("C:/BiS335-Final-Project/Dataset/midterm/mutation.rds")

gex <- na.omit(gex)
clin <- na.omit(clin)

id <- clin$sample_id

exp_gene <- rownames(gex) 
exp_sample <- colnames(gex) 

with_mut_id <- mut$sample_id # sample ids that have a mutation 
unique_mut_gene <- unique(mut$Hugo_Symbol) # Remove redundant gene names by unique() function. 
mut_hugo <- mut$Hugo_Symbol

common_pool <- intersect(intersect(id, with_mut_id), exp_sample)

stage <- clin$stage[match(common_pool, id)]
subtype <- clin$subtype[match(common_pool, id)]
surv <- clin$survival_time[match(common_pool, id)]
gex <- gex[,match(common_pool, exp_sample)]

####################
## 1. Expression 
####################

getSigmoidGene <- function(surv, gex, a) {
  
  # a : moving average (smoothing) width
  r <- c()
  rs <-c() # save correlation coeffiecient
  ps <-c() # save p-value of correlation test
  score <- slope <- c() # save score and slope of legit
  
    for (i in 1:nrow(gex)){
    gene_expr <- gex[i,] # i th gene expression 
    
    temp <- surv[length(surv)+1-order(gene_expr, decreasing = TRUE)] 
    temp <- convolve(temp, matrix(1/a, ncol = a), type = "filter")
    temp2 <- convolve(sort(gene_expr), matrix(1/a, ncol = a), type = "filter")
    
    norm <- (temp-min(temp)+0.1) # scaling
    norm <- norm/(max(norm)+1)
    
    legit <- log(norm/(1-norm))
    
    r[i] <- cor(legit, temp2)
    slope[i] <- (mean(legit[1:round(length(legit)/2)])-mean(legit[round(length(legit)/2):length(legit)]))/
      (mean(temp2[1:round(length(legit)/2)])-mean(temp2[round(length(legit)/2):length(legit)]))
    score[i] <- abs(slope[i]*r[i])
    
    t <- r[i] * sqrt((length(legit) - 2)/(1 - r[i]^2))
    p <- min(1, 2 * pt(t, length(legit) - 2))
    if(is.na(p)){ #don`t save when p is NA
      ps[i] <- 1
    }
    else{
      ps[i] <- p
    }
    if(!is.na(r[i])){ #don`t save when r is NA
      rs[i] <- r[i]
    }
  }
  
  # match(sort(r[r>0], decreasing = TRUE), r)
  # match(sort(r[r<0]), r)
  
  data <- data.frame(cor = r[order(score, decreasing = TRUE)])
  rownames(data) = exp_gene[order(score, decreasing = TRUE)]
  data <- cbind(data, slope = slope[order(score, decreasing = TRUE)])
  data <- cbind(data, number = c(1:length(r))[order(score, decreasing = TRUE)])
  data <- cbind(data, score = score[order(score, decreasing = TRUE)])
  
  return(data)
}

sig.g.1st <- getSigmoidGene(surv, gex, 5)
write.csv(sig.g.1st, file = "sigmoid_1_a5.csv", row.names = T)

plotSigmoidGene <- function(surv, gex, a, data) {
  
  # a : moving average (smoothing) width
  r <- c()
  rs <-c() #save correlation coeffiecient
  ps <-c() #save p-value of correlation test
  score <- slope <- c() # save score and slope of legit
  
  par(mfrow = c(3, 5))
  for (i in 1:15) {
    gene_expr <- gex[data[i, 3],] # i th gene expression
    
    temp <- surv[length(surv)+1-order(gene_expr, decreasing = TRUE)]
    temp <- convolve(temp, matrix(1/a, ncol = a), type = "filter")
    temp2 <- convolve(sort(gene_expr), matrix(1/a, ncol = a), type = "filter")
    
    norm <- (temp-min(temp)+0.1) # scaling
    norm <- norm/(max(norm)+1)
    
    legit <- log(norm/(1-norm))
    
    r[i] <- cor(legit, temp2)
    slope[i] <- (mean(legit[1:round(length(legit)/2)])-mean(legit[round(length(legit)/2):length(legit)]))/
      (mean(temp2[1:round(length(legit)/2)])-mean(temp2[round(length(legit)/2):length(legit)]))
    score[i] <- abs(slope[i]*r[i])
    
    t <- r[i] * sqrt((length(legit) - 2)/(1 - r[i]^2))
    p <- min(1, 2 * pt(t, length(legit) - 2))
    if(is.na(p)){ #don`t save when p is NA
      ps[i] <- 1
    }
    else{
      ps[i] <- p
    }
    if(!is.na(r[i])){ #don`t save when r is NA
      rs[i] <- r[i]
    }
    plot(temp2, legit, cex = 0.5, main = rownames(data)[i], xlab = "gene exp. level")
    text(x = min(temp2), y = 5, labels = paste0("score = ", round(data[i, 4], 2), "(", round(data[i, 1], 2), ", ", round(data[i, 2], 2), ")"), pos = 4, col = 'blue')
  }
}

plotSigmoidGene(surv, gex, 5, sig.g.1st)
# 
# # Plot gexp. - legit
# 
# a = 5 # moving average (smoothing) width
# r <- c()
# rs <-c() #save correlation coeffiecient
# ps <-c() #save p-value of correlation test
# score <- slope <- c() # save score and slope of legit
# 
# par(mfrow = c(3, 5))
# for (i in 1:15) {
#   gene_expr <- gex[sig.g.1st$number[i],] # i th gene expression
# 
#   temp <- surv[length(surv)+1-order(gene_expr, decreasing = TRUE)]
#   temp <- convolve(temp, matrix(1/a, ncol = a), type = "filter")
#   temp2 <- convolve(sort(gene_expr), matrix(1/a, ncol = a), type = "filter")
# 
#   norm <- (temp-min(temp)+0.1) # scaling
#   norm <- norm/(max(norm)+1)
# 
#   legit <- log(norm/(1-norm))
# 
#   r[i] <- cor(legit, temp2)
#   slope[i] <- (mean(legit[1:round(length(legit)/2)])-mean(legit[round(length(legit)/2):length(legit)]))/
#     (mean(temp2[1:round(length(legit)/2)])-mean(temp2[round(length(legit)/2):length(legit)]))
#   score[i] <- abs(slope[i]*r[i])
# 
#   t <- r[i] * sqrt((length(legit) - 2)/(1 - r[i]^2))
#   p <- min(1, 2 * pt(t, length(legit) - 2))
#   if(is.na(p)){ #don`t save when p is NA
#     ps[i] <- 1
#   }
#   else{
#     ps[i] <- p
#   }
#   if(!is.na(r[i])){ #don`t save when r is NA
#     rs[i] <- r[i]
#   }
#   plot(temp2, legit, cex = 0.5, main = rownames(sig.g.1st)[i], xlab = "gene exp. level")
#   text(x = min(temp2), y = 5, labels = paste0("score = ", round(sig.g.1st$score[i], 2), "(", round(sig.g.1st$cor[i], 2), ", ", round(sig.g.1st$slope[i], 2), ")"), pos = 4, col = 'blue')
# }


# sum((surv-mean(surv))^2)/length(surv)
# norm <- (surv+0.1)/max(surv+1.1)
# legit <- log(norm/(1-norm))

# psm <- ps[ps<0.05]
# #adjust the data with multiplication correction method by "holm`s method"`
# ps_h <- p.adjust(ps, method = "bonferroni")
# psm_holm <- ps[ps_h<0.05]
# gene_cor <- rownames(gex)[match(psm, ps)] #genes having p value lower then 0.05
# gene_cor_h <- rownames(gex)[match(psm_holm, ps)] #genes having p value lower then 0.05 even after holm method


getMSEGene <- function(surv, gex, idx) {
  MSE <- c()
  minMSE <- c()
  thr_idx <- c()
  thr_exp <- c()
  
  par(mfrow = c(3, 5), ps = 15)
  
  for (j in idx) {
    MSE <- c()
    value <- sort(unique(gex[j,]))
    value <- convolve(value, c(1/2,1/2), type = "filter")
    # sizelim <- floor(length(surv)*0.05)
    for (i in 20:(length(value)-20)){
      group1 <- surv[gex[j,]>value[i]]
      group2 <- surv[gex[j,]<value[i]]
      MSE[i-19] <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
    }
    plot(MSE, cex = 0.5, main = paste0(rownames(gex)[j], " (", which.min(MSE)+19, ", ", round(min(MSE), 0), ")"))
    abline(h = sum((surv-mean(surv))^2)/length(surv), col = "green", lty = 2)
    abline(h = min(MSE), col = "red", lty = 2)
    abline(v = which.min(MSE), col="red", lty = 1)
    minMSE <- c(minMSE, min(MSE))
    thr_idx <- c(thr_idx, which.min(MSE)+19)
    thr_exp <- c(thr_exp, value[which.min(MSE)+19])
  }
  
  data <- data.frame(idx = thr_idx, exp = thr_exp, MSE = minMSE)
  rownames(data) = rownames(gex)[idx]

  return(data)
}

MSE.g.1st <- getMSEGene(surv, gex, sig.g.1st$number[1:15])
write.csv(MSE.g.1st, file = "MSE_gene_1.csv", row.names = TRUE)



####################
## 2. Mutation 
####################

getMSEMut <- function(surv, gene) {
  MSE <- c()
  O <- c()
  X <- c()
  
  for (i in 1:length(unique_mut_gene)) {
    mutO <- with_mut_id[mut_hugo == unique_mut_gene[i]] # id of samples with i th gene mutation  
    mutX <- setdiff(with_mut_id, mutO) # id of samples without i th gene mutation
    group1 <- na.omit(surv[match(mutO, common_pool)])
    group2 <- na.omit(surv[match(mutX, common_pool)])
    MSE <- c(MSE, (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2)))
    O <- c(O, length(mutO))
    X <- c(X, length(mutX))
  }
  
  data <- data.frame(mutO = O, mutX = X, MSE = MSE) 
  rownames(data) = unique_mut_gene
  data <- data[order(MSE), ]
  
  return(data)
}

MSE.m.1st <- getMSEMut(surv, unique_mut_gene)
MSE.m.1st.filter <- MSE.m.1st[MSE.m.1st$mutO>=20,]
write.csv(MSE.m.1st.filter, file = "MSE_mut_1_filtered.csv", row.names = TRUE)

plotMSEMut <- function(surv, gene, mut.id, MSE.data) {
  par(mfrow = c(3, 5))
  for(i in 1:15) {
    mutO <- mut.id[mut_hugo == rownames(MSE.data)[i]]
    mutX <- setdiff(mut.id, mutO) 
    group1 <- na.omit(surv[match(mutO, common_pool)])
    group2 <- na.omit(surv[match(mutX, common_pool)])
    MSE <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
    data <- data.frame(surv = c(group1, group2), group = c(rep("mutO", length(group1)), rep("mutX", length(group2))))
    boxplot(surv ~ group, data = data, main = paste0(rownames(MSE.data)[i], "(", length(mutO), ", ", round(MSE, 0), ")"))  
  }
}

plotMSEMut(surv, gene, with_mut_id, MSE.m.1st.filter)

# # Boxplot mut - surv
# par(mfrow = c(3, 5))
# for(i in 1:15) {
#   mutO <- with_mut_id[mut_hugo == rownames(MSE.m.1st.filter)[i]]
#   mutX <- setdiff(with_mut_id, mutO) 
#   group1 <- na.omit(surv[match(mutO, common_pool)])
#   group2 <- na.omit(surv[match(mutX, common_pool)])
#   MSE <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
#   data <- data.frame(surv = c(group1, group2), group = c(rep("mutO", length(group1)), rep("mutX", length(group2))))
#   boxplot(surv ~ group, data = data, main = paste0(rownames(MSE.m.1st.filter)[i], "(", length(mutO), ", ", round(MSE, 0), ")"))  
# }