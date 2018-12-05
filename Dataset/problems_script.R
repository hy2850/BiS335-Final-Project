rm(list = ls()) # Clear all workspace
setwd("C:/Users/박현철/Desktop/BiS335 바이오통계/Final Project/Dataset")

## Environment settings
# setwd to directory containing below rds files!
clin <- readRDS("clinical.rds")
mut <- readRDS("mutation.rds")
gex <- readRDS("expression.rds")

library(dplyr)

#############################################################################
library(survival)
# library(survminer)

#subtype <- clin$subtype


death_time = sample(1:50, 100, replace=TRUE)
status = sample(1:0, 100, replace=TRUE)
group = c(rep(0,50), rep(1, 50))

group <- factor(group)                        # 0, 1 2개의 level로 factor화
fit = survfit(Surv(death_time, status)~group)

plot( fit, ylab="Survival", xlab="Days", col=c("red","blue"), lty=1:2, mark.time=T)
legend("topright", legend=levels(group), col=c("red","blue"), lty=1:2)
# 출처: http://sosal.kr/865 [so_sal　]

data(colon)
attach(colon)
str(colon)