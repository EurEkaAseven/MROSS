rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(MASS)
library(kerndwd)
source("basic/DataGenerator.R")
source("basic/getMLE.R")

nrep <-200                  # Replication size
N <- 5e5                    # Full data size
p <- 7                      # Dimension of X
r0 <- 1000                  # Step 1 sampling size
lambda = 0
kern = vanilladot()
# Generating model
gmodel <- "logistic"
option <- "i"            # Distribution setting of X
theta_true <- c(0.5,rep(0.5,p))
DATA <- GenerateData(N = 10*N, p = p, theta = theta_true, 
                     option = option, model = gmodel)
theta_true_dwd <- as.vector(kerndwd(DATA[,-c(1,2)], DATA[,1], 
                            kern = kern,
                            qval = 1, lambda = lambda,
                            eps = 1e-5, maxit = 1e+5)$alpha)
rm("DATA")


cputime_log <- vector()
cputime_dwd <- vector()
theta_all_log <- vector()
theta_all_dwd <- vector()
Time_rec <- matrix(nrow = 2, ncol = 4)
Time_rec <- as.data.frame(Time_rec)
colnames(Time_rec) <- c("mean", "median", "sd", "max")
rownames(Time_rec) <- c("logistic8","dwd8")

EMSE <- Time_rec[,1:2]
colnames(EMSE) <- c("MSE","sdMSE")
for (i in 1:nrep) {
  cat(i)
  set.seed(30231102+i)
  DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
  gc()
  init.time <- Sys.time()
  theta_new <- getMLE(DATA[,-1], DATA[,1])$par
  end.time <- Sys.time()
  cputime <- difftime(end.time, init.time, units = "secs")
  
  cputime_log <- c(cputime_log, cputime)
  theta_all_log <- cbind(theta_all_log, theta_new)
  
  gc()
  init.time <- Sys.time()
  theta_new <- kerndwd(DATA[,-c(1,2)], DATA[,1], 
                       kern = kern,
                       qval = 1, lambda = lambda,
                       eps = 1e-5, maxit = 1e+5)$alpha
  end.time <- Sys.time()
  cputime <- difftime(end.time, init.time, units = "secs")
  cputime_dwd <- c(cputime_dwd, cputime)
  theta_all_dwd <- cbind(theta_all_dwd, c(theta_new))
  rm("DATA")
}
Time_rec[1,] <- c(mean(cputime_log), median(cputime_log), sd(cputime_log), max(cputime_log))
Time_rec[2,] <- c(mean(cputime_dwd), median(cputime_dwd), sd(cputime_dwd), max(cputime_dwd))

  bias <- (theta_all_log-theta_true)^2
  MSE <- colSums(bias)
  EMSE[1,1] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  EMSE[1,2] <- sd(MSE, na.rm = T)
  
  bias <- (theta_all_dwd-theta_true_dwd)^2
  MSE <- colSums(bias)
  EMSE[2,1] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  EMSE[2,2] <- sd(MSE, na.rm = T)
  
  save(list = ls(), file=paste0("FullTime_p=",p,".RData"))


  
rm(list = ls())
gc()
source("DataGenerator.R")
source("getMLE.R")

nrep <-200                   # Replication size
N <- 5e5                     # Full data size
p <- 20                      # Dimension of X
r0 <- 1000                   # Step 1 sampling size
lambda = 0
kern = vanilladot()
# Generating model
gmodel <- "logistic"
option <- "i"            # Distribution setting of X
theta_true <- c(0.5,rep(0.5,p))
DATA <- GenerateData(N = 10*N, p = p, theta = theta_true, 
                     option = option, model = gmodel)
theta_true_dwd <- as.vector(kerndwd(DATA[,-c(1,2)], DATA[,1], 
                                    kern = kern,
                                    qval = 1, lambda = lambda,
                                    eps = 1e-5, maxit = 1e+5)$alpha)
rm("DATA")


cputime_log <- vector()
cputime_dwd <- vector()
theta_all_log <- vector()
theta_all_dwd <- vector()
Time_rec <- matrix(nrow = 2, ncol = 4)
Time_rec <- as.data.frame(Time_rec)
colnames(Time_rec) <- c("mean", "median", "sd", "max")
rownames(Time_rec) <- c("logistic21","dwd21")

EMSE <- Time_rec[,1:2]
colnames(EMSE) <- c("MSE","sdMSE")
for (i in 1:nrep) {
  cat(i)
  set.seed(30231102+i)
  DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
  gc()
  init.time <- Sys.time()
  theta_new <- getMLE(DATA[,-1], DATA[,1])$par
  end.time <- Sys.time()
  cputime <- difftime(end.time, init.time, units = "secs")
  
  cputime_log <- c(cputime_log, cputime)
  theta_all_log <- cbind(theta_all_log, theta_new)
  
  gc()
  init.time <- Sys.time()
  theta_new <- kerndwd(DATA[,-c(1,2)], DATA[,1], 
                       kern = kern,
                       qval = 1, lambda = lambda,
                       eps = 1e-5, maxit = 1e+5)$alpha
  end.time <- Sys.time()
  cputime <- difftime(end.time, init.time, units = "secs")
  cputime_dwd <- c(cputime_dwd, cputime)
  theta_all_dwd <- cbind(theta_all_dwd, c(theta_new))
  rm("DATA")
}
Time_rec[1,] <- c(mean(cputime_log), median(cputime_log), sd(cputime_log), max(cputime_log))
Time_rec[2,] <- c(mean(cputime_dwd), median(cputime_dwd), sd(cputime_dwd), max(cputime_dwd))

bias <- (theta_all_log-theta_true)^2
MSE <- colSums(bias)
EMSE[1,1] <- mean(MSE, na.rm = T)
cat(sum(is.na(MSE)))
cat(" ")
EMSE[1,2] <- sd(MSE, na.rm = T)

bias <- (theta_all_dwd-theta_true_dwd)^2
MSE <- colSums(bias)
EMSE[2,1] <- mean(MSE, na.rm = T)
cat(sum(is.na(MSE)))
cat(" ")
EMSE[2,2] <- sd(MSE, na.rm = T)

save(list = ls(), file=paste0("FullTime_p=",p,".RData"))
