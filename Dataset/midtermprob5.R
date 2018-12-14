# put the appropriate path name in your computer
clin <- readRDS("C:/Users/Lee/Desktop/Midterm-Project/Midterm-Project/clin.rds")
mut <- readRDS("C:/Users/Lee/Desktop/Midterm-Project/Midterm-Project/mutation.rds")
gex <- readRDS("C:/Users/Lee/Desktop/Midterm-Project/Midterm-Project/expression.rds")

gex <- na.omit(gex)
clin <- na.omit(clin)

id <- clin$sample_id
with_mut_id <- mut$sample_id        # sample ids that have a mutation 
unique_mut_gene <- unique(mut$Hugo_Symbol) # Remove redundant gene names by unique() function. 
mut_hugo <- mut$Hugo_Symbol 

exp_gene <- rownames(gex) 
exp_sample <- colnames(gex) 

common_pool <- intersect(intersect(id, with_mut_id), exp_sample)

stage <- clin$stage[match(common_pool, id)]
subtype <- clin$subtype[match(common_pool, id)]
surv <- clin$survival_time[match(common_pool, id)]
gex <- gex[,match(common_pool, exp_sample)]

a = 10 # moving average (smoothing) width 
r <- c()
rs <-c() #save correlation coeffiecient
ps <-c() #save p-value of correlation test
r <- c()
for (i in 1:nrow(gex)){
  gene_expr <- gex[i,] # i th gene expression 
  
  temp <- surv[length(surv)+1-order(gene_expr, decreasing = TRUE)] 
  temp <- convolve(temp, matrix(1/a, ncol = a), type = "filter")
  
  norm <- (temp-min(temp)+0.1) # scaling
  norm <- norm/(max(norm)+1)
  
  legit <- log(norm/(1-norm))
  
  r[i] <- cor(legit, c(1:length(legit))/length(legit))
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
  # plot(legit, cex = 0.5, main = i)
  # Sys.sleep(1)
}

match(sort(r[r>0], decreasing = TRUE), r)
match(sort(r[r<0]), r)

data <- data.frame(cor = r[order(abs(r), decreasing = TRUE)])
rownames(data) = exp_gene[order(abs(r), decreasing = TRUE)]
data <- cbind(data, number = c(1:length(r))[order(abs(r), decreasing = TRUE)])

data

# 
sum((surv-mean(surv))^2)/length(surv)
# norm <- (surv+0.1)/max(surv+1.1)
# legit <- log(norm/(1-norm))
# 
# 
# 
# rs <-c() #save correlation coeffiecient
# ps <-c() #save p-value of correlation test
# r <- c()
# for (i in 1: length(rownames(gex))){
#   gene_expr <- gex[i,]
#   r[i] <- cor(gene_expr, legit)
#   t <- r[i] * sqrt((length(legit) - 2)/(1 - r^2))
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
# }
# 
psm <- ps[ps<0.05]
#adjust the data with multiplication correction method by "holm`s method"`
ps_h <- p.adjust(ps, method = "bonferroni")
psm_holm <- ps[ps_h<0.05]
gene_cor <- rownames(gex)[match(psm, ps)] #genes having p value lower then 0.05
gene_cor_h <- rownames(gex)[match(psm_holm, ps)] #genes having p value lower then 0.05 even after holm method

j = 10656
MSE <- c()
value <- sort(unique(gex[j,]))
value <- convolve(value, c(1/2,1/2), type = "filter")
for (i in 20:(length(value)-20)){
  group1 <- surv[gex[j,]>value[i]]
  group2 <- surv[gex[j,]<value[i]]
  MSE[i-19] <- sum((group1-mean(group1))^2)/length(group1)+sum((group2-mean(group2))^2)/length(group2)
  MSE[i-19] <- MSE[i-19]/2
}
plot(MSE, cex = 0.5)
min(MSE)
