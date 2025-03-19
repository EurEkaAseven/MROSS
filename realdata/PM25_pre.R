rm(list = ls())
gc()
######## set your working directory here ########
setwd("/Users/yujun/Desktop/realdata_code_4.21")
#################################################
library(kerndwd)
library(matrixStats)
library(data.table)
source("getMLE.R")

# read the data
# PM2.5 data
data <- NULL
# change the path to where you store the data
for(i in list.files("/Users/yujun/Documents/PM25")){
  data1 <- read.csv(paste0("/Users/yujun/Documents/PM25/",i))
  data1 <- as.matrix(data1[,c(6,8:14,17)])
  data1 <- cbind(c(data1[-1,1],NA),data1[,2:9]) 
  data <- rbind(data,data1)
}
data <- data[complete.cases(data),]


Y <- ifelse(data[,1]>75,1,0)#/sd(data$PM2.5)
X <- as.matrix((data[,-1]))
X <- cbind(1,X)#scale(model.matrix(A)[,-1],center = T)




DATA <- data.table(Y = Y,
                   X = X)
DATA <- as.matrix(DATA)
N <- nrow(DATA)

ntest = ceiling(0.2*nrow(DATA))
set.seed(0)
sample_index <- sample(N, ntest, replace = FALSE)
testdata <- DATA[sample_index,]
DATA <- DATA[-sample_index,]


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
save(DATA, testdata, N, p, theta_dwd, theta_logistic, file="PM25_pre.Rdata")
