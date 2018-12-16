rm(list = ls()) # Clear all workspace


## Environment settings
# setwd to directory containing below rds files!
clin <- readRDS("clinical.rds")
mut <- readRDS("mutation.rds")
gex <- readRDS("expression.rds")

library(survival)
library(survminer)

# subtype <- clin$subtype
# sam_id <- clin$sample_id[!is.na(subtype)]  
# surv_t <- clin$survival_time[!is.na(subtype)]  
# surv_ind <- clin$survival_index[!is.na(subtype)]
# vital <- clin$vital_status[!is.na(subtype)]
# stage <- clin$stage[!is.na(subtype)] 
# subtype <- subtype[!is.na(subtype)]


clin <- na.omit(clin) # remove NA including data
surv_t <- clin$survival_time
vital <- clin$vital_status
stage <- clin$stage
subtype <- clin$subtype
surv_ind <- clin$survival_index
vital <- as.numeric(vital)

su <- Surv(surv_t, vital) # test for stage
s_diff <- survdiff(su ~ stage)
s_diff

s_diff <- survdiff(su ~ subtype)# test for subtype
s_diff

s_diff <- survdiff(su ~ surv_ind)# test for surv index
s_diff

fit0 <- survfit(su~surv_ind, data = clin) # plotting survival curve
ggsurvplot(fit0)

#-------------------common sample sets ------------------------

common_id <- intersect(intersect(clin$sample_id, mut$sample_id), colnames(gex))
data2 <- clin[which(clin$sample_id %in% common_id), 1:6]

su2 <- Surv(data2$survival_time, as.numeric(data2$vital_status)) # test for stage
s_diff <- survdiff(su2 ~ data2$stage)
s_diff

s_diff <- survdiff(su2 ~ data2$subtype)# test for subtype
s_diff

s_diff <- survdiff(su2 ~ data2$survival_index)# test for surv index
s_diff

#------------------1.d------------------

surv1 <- surv_t[surv_ind == 1]
surv2 <- surv_t[surv_ind == 2]
surv3 <- surv_t[surv_ind == 3]
surv4 <- surv_t[surv_ind == 4]
stage1 <- surv_t[stage == sort(unique(stage))[1]]
stage2 <- surv_t[stage == sort(unique(stage))[2]]
stage3 <- surv_t[stage == sort(unique(stage))[3]]

median(surv1); median(surv2); median(surv3); median(surv4)
median(stage1); median(stage2); median(stage3)

par(mfrow = c(1,2))
barplot(c(median(surv1), median(surv2), median(surv3), median(surv4)), 
        xlab = "survival index",
        ylab = "survival time",
        main = "median for each survival index", names.arg = c(1, 2, 3, 4))
barplot(c(median(stage1), median(stage2), median(stage3)), 
        xlab = "stage",
        ylab = "survival time",
        main = "median for each stage", names.arg = c(1, 2, 3))
