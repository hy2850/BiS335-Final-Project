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

# sum((surv-mean(surv))^2)/length(surv)
# norm <- (surv+0.1)/max(surv+1.1)
# legit <- log(norm/(1-norm))

# psm <- ps[ps<0.05]
# #adjust the data with multiplication correction method by "holm`s method"`
# ps_h <- p.adjust(ps, method = "bonferroni")
# psm_holm <- ps[ps_h<0.05]
# gene_cor <- rownames(gex)[match(psm, ps)] #genes having p value lower then 0.05
# gene_cor_h <- rownames(gex)[match(psm_holm, ps)] #genes having p value lower then 0.05 even after holm method

getMSEGene <- function(surv, gex) {
  MSE <- c()
  minMSE <- c()
  thr_idx <- c()
  thr_exp <- c()
  
  for (j in 1:dim(gex)[1]) {
    MSE <- c()
    value <- sort(unique(gex[j,]))
    value <- convolve(value, c(1/2,1/2), type = "filter")
    # sizelim <- floor(length(surv)*0.05)
    for (i in 20:(length(value)-20)){
      group1 <- surv[gex[j,]>value[i]]
      group2 <- surv[gex[j,]<value[i]]
      MSE[i-19] <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
    }
    minMSE <- c(minMSE, min(MSE))
    thr_idx <- c(thr_idx, which.min(MSE)+19)
    thr_exp <- c(thr_exp, value[which.min(MSE)+19])
  }
  
  data <- data.frame(number = c(1:dim(gex)[1]), idx = thr_idx, exp = thr_exp, MSE = minMSE)
  rownames(data) = rownames(gex)
  
  return(data)
}

plotMSEGene <- function(surv, gex, data) {
  MSE <- c()

  par(mfrow = c(3, 3), ps = 15)

  for (j in 1:9) {
    MSE <- c()
    value <- sort(unique(gex[data[j, 1],]))
    value <- convolve(value, c(1/2,1/2), type = "filter")

    for (i in 20:(length(value)-20)){
      group1 <- surv[gex[data[j, 1],]>value[i]]
      group2 <- surv[gex[data[j, 1],]<value[i]]
      MSE[i-19] <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
    }
    plot(MSE, cex = 0.5, main = paste0(rownames(gex)[data[j, 1]], " (", which.min(MSE)+19, ", ", round(min(MSE), 0), ")"))
    abline(h = sum((surv-mean(surv))^2)/length(surv), col = "green", lty = 2)
    abline(h = min(MSE), col = "red", lty = 2)
    abline(v = which.min(MSE), col="red", lty = 1)
  }
}



####################
## 2. Mutation 
####################

getMSEMut <- function(surv, pool) {
  MSE <- c()
  O <- c()
  X <- c()
  
  for (i in 1:length(unique_mut_gene)) {
    mutO <- intersect(with_mut_id[mut_hugo == unique_mut_gene[i]], pool) # id of samples with i th gene mutation  
    mutX <- setdiff(pool, mutO) # id of samples without i th gene mutation
    group1 <- na.omit(surv[match(mutO, pool)])
    group2 <- na.omit(surv[match(mutX, pool)])
    MSE <- c(MSE, (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2)))
    O <- c(O, length(mutO))
    X <- c(X, length(mutX))
  }
  
  data <- data.frame(mutO = O, mutX = X, MSE = MSE) 
  rownames(data) = unique_mut_gene
  data <- data[order(MSE), ]
  
  return(data)
}

plotMSEMut <- function(surv, pool, MSE.data) {
  par(mfrow = c(3, 3))
  for(i in 1:9) {
    mutO <- intersect(with_mut_id[mut_hugo == rownames(MSE.data)[i]], pool)
    mutX <- setdiff(pool, mutO) 
    group1 <- surv[match(mutO, pool)]
    group2 <- surv[match(mutX, pool)]
    MSE <- (sum((group1-mean(group1))^2)+sum((group2-mean(group2))^2))/(length(group1)+length(group2))
    data <- data.frame(surv = c(group1, group2), group = c(rep("mutO", length(group1)), rep("mutX", length(group2))))
    boxplot(surv ~ group, data = data, main = paste0(rownames(MSE.data)[i], "(", length(mutO), ", ", round(MSE, 0), ")"))  
  }
}



####################
## 3. Grouping
####################

# 1st Grouping
# MSE before grouping = 1953133

# for gene expression data
MSE.g.1st <- getMSEGene(surv, gex)
MSE.g.1st <- MSE.g.1st[order(MSE.g.1st$MSE),]
write.csv(MSE.g.1st, file = "MSE_gene_1.csv", row.names = TRUE)

plotMSEGene(surv, gex, MSE.g.1st)
# 1st minMSE = 1783938 (ZFHX3, -0.8195, 26) 
# 2nd minMSE = 1788354 (PELI1, 0.57325, 218)

# for mutation data
MSE.m.1st <- getMSEMut(surv, common_pool)
#MSE.m.1st.filter <- MSE.m.1st[MSE.m.1st$mutO>=20,]
write.csv(MSE.m.1st, file = "MSE_mut_1.csv", row.names = TRUE)

plotMSEMut(surv, common_pool, MSE.m.1st)
# minMSE = 1840273 (TTN, 360, 784)

# using ZFHX3 -0.8195 as a threshold
ZFHX3.idx <- MSE.g.1st[1, 1]
ZFHX3.thr <- MSE.g.1st[1, 3]
L.surv <- surv[gex[ZFHX3.idx,]<ZFHX3.thr] # 26 obs.
R.surv <- surv[gex[ZFHX3.idx,]>ZFHX3.thr] # 404 obs.
L.gex <- gex[, gex[ZFHX3.idx,]<ZFHX3.thr]
R.gex <- gex[, gex[ZFHX3.idx,]>ZFHX3.thr]
L.pool <- colnames(gex)[gex[ZFHX3.idx,]<ZFHX3.thr]
R.pool <- colnames(gex)[gex[ZFHX3.idx,]>ZFHX3.thr]


# 2nd Grouping (for ZFHX3 R group)
# MSE before grouping = 1455782

# for gene expression data
MSE.g.2nd <- getMSEGene(R.surv, R.gex)
MSE.g.2nd <- MSE.g.2nd[order(MSE.g.2nd$MSE),]
write.csv(MSE.g.2nd, file = "MSE_gene_2.csv", row.names = TRUE)

plotMSEGene(R.surv, R.gex, MSE.g.2nd)
# 1st minMSE = 1328730 (ZFP36, 0.92575, 142) 

# for mutation data
MSE.m.2nd <- getMSEMut(R.surv, R.pool)
# MSE.m.2nd.filter <- MSE.m.2nd[MSE.m.2nd$mutO>=20,]
write.csv(MSE.m.2nd, file = "MSE_mut_2.csv", row.names = TRUE)

plotMSEMut(R.surv, R.pool, MSE.m.2nd)
# minMSE = 1322645 (MPP7, 1, 403)

# using ZFP36 0.92575 as a threshold
ZFP36.idx <- MSE.g.2nd[1, 1]
ZFP36.thr <- MSE.g.2nd[1, 3]
RL.surv <- R.surv[R.gex[ZFP36.idx,]<ZFP36.thr] # 146 obs.
RR.surv <- R.surv[R.gex[ZFP36.idx,]>ZFP36.thr] # 258 obs.
RL.gex <- R.gex[, R.gex[ZFP36.idx,]<ZFP36.thr]
RR.gex <- R.gex[, R.gex[ZFP36.idx,]>ZFP36.thr]
RL.pool <- colnames(R.gex)[R.gex[ZFP36.idx,]<ZFP36.thr]
RR.pool <- colnames(R.gex)[R.gex[ZFP36.idx,]>ZFP36.thr]


# 3rd Grouping (1) (for ZFP36 L group)
# MSE before grouping = 738770.3

# for gene expression data
MSE.g.3rd.1 <- getMSEGene(RL.surv, RL.gex)
MSE.g.3rd.1 <- MSE.g.3rd.1[order(MSE.g.3rd.1$MSE),]
write.csv(MSE.g.3rd.1, file = "MSE_gene_3_1.csv", row.names = TRUE)

plotMSEGene(RL.surv, RL.gex, MSE.g.3rd.1)
# minMSE = 589888 (FCGR1A, 3.09775, 31) 

# for mutation data
MSE.m.3rd.1 <- getMSEMut(RL.surv, RL.pool)
write.csv(MSE.m.3rd.1, file = "MSE_mut_3_1.csv", row.names = TRUE)

plotMSEMut(RL.surv, RL.pool, MSE.m.3rd.1)
# minMSE = 567615 (OR5W2, 1, 145)

# using FCGR1A 3.09775 as a threshold
FCGR1A.idx <- MSE.g.3rd.1[1, 1]
FCGR1A.thr <- MSE.g.3rd.1[1, 3]
RLL.surv <- RL.surv[RL.gex[FCGR1A.idx,]<FCGR1A.thr] # 32 obs.
RLR.surv <- RL.surv[RL.gex[FCGR1A.idx,]>FCGR1A.thr] # 114 obs.
RLL.gex <- RL.gex[, RL.gex[FCGR1A.idx,]<FCGR1A.thr]
RLR.gex <- RL.gex[, RL.gex[FCGR1A.idx,]>FCGR1A.thr]
RLL.pool <- colnames(RL.gex)[RL.gex[FCGR1A.idx,]<FCGR1A.thr]
RLR.pool <- colnames(RL.gex)[RL.gex[FCGR1A.idx,]>FCGR1A.thr]


# 3rd Grouping (2) (for ZFP36 R group)
# MSE before grouping = 1662584

# for gene expression data
MSE.g.3rd.2 <- getMSEGene(RR.surv, RR.gex)
MSE.g.3rd.2 <- MSE.g.3rd.2[order(MSE.g.3rd.2$MSE),]
write.csv(MSE.g.3rd.2, file = "MSE_gene_3_2.csv", row.names = TRUE)

plotMSEGene(RR.surv, RR.gex, MSE.g.3rd.2)
# minMSE = 1480910 (DHRS9, -1.069125, 25) 

# for mutation data
MSE.m.3rd.2 <- getMSEMut(RR.surv, RR.pool)
write.csv(MSE.m.3rd.2, file = "MSE_mut_3_2.csv", row.names = TRUE)

plotMSEMut(RR.surv, RR.pool, MSE.m.3rd.2)
# minMSE = 567615 (MPP7, 1, 257)

# using FCGR1A 3.09775 as a threshold
FCGR1A.idx <- MSE.g.3rd.1[1, 1]
FCGR1A.thr <- MSE.g.3rd.1[1, 3]
RLL.surv <- RL.surv[RL.gex[FCGR1A.idx,]<FCGR1A.thr] # 32 obs.
RLR.surv <- RL.surv[RL.gex[FCGR1A.idx,]>FCGR1A.thr] # 114 obs.
RLL.gex <- RL.gex[, RL.gex[FCGR1A.idx,]<FCGR1A.thr]
RLR.gex <- RL.gex[, RL.gex[FCGR1A.idx,]>FCGR1A.thr]
RLL.pool <- colnames(RL.gex)[RL.gex[FCGR1A.idx,]<FCGR1A.thr]
RLR.pool <- colnames(RL.gex)[RL.gex[FCGR1A.idx,]>FCGR1A.thr]

