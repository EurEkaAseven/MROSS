rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility/realdata")
#################################################
library(kerndwd)
library(matrixStats)
library(data.table)
source("getMLE.R")

# read the data
DATA <- fread("SUSY.csv", header = FALSE)
DATA <- data.table(Y = DATA[,1],
                   X = cbind(1,DATA[,-1]))
DATA <- as.matrix(DATA)

testdata <- tail(DATA, 5e5)
DATA <- head(DATA,nrow(DATA)-5e5)

N <- nrow(DATA)
p <- ncol(DATA)-2

# scaling
means <- colMeans(DATA)
sds <- colSds(DATA)

DATA[,-c(1,2)] <- t((t(DATA[,-c(1,2)]) - means[-c(1,2)])/sds[-c(1,2)])
testdata[,-c(1,2)] <- t((t(testdata[,-c(1,2)]) - means[-c(1,2)])/sds[-c(1,2)])

# theta_full
theta_dwd <- as.vector(kerndwd(DATA[,-c(1,2)], DATA[,1], 
                               kern = vanilladot(),
                               qval = 1, lambda = 0,
                               eps = 1e-5, maxit = 1e+5)$alpha)
theta_logistic <- as.vector(getMLE(DATA[,-1], DATA[,1])$par)

# saving
save(DATA, testdata, N, p, theta_dwd, theta_logistic, file="SUSY_pre.Rdata")
