rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Reply1")
#################################################
library(MASS)
library(Rcpp)
source("basic/getMLE.R")
sourceCpp("Rcpp/OSMAC.cpp")
sourceCpp("Rcpp/RBIMP.cpp")

load("realdata/SUSY_pre.Rdata")
rm("testdata")
gc()

nrep <-200                 # Replication size
r0 <- 1000                      # Step 1 sampling size
r.ss <- c(2000, 5000)
# Working model
wmodel <- "logistic"

EMSE <- matrix(NA,nrow = length(r.ss),ncol = 2)
EMSE <- as.data.frame(EMSE)
rownames(EMSE) <- r.ss
colnames(EMSE) <- c("OSMAC","RBIMP")
sdMSE <- EMSE

############## OSMAC ################
#r.ss <- c(2000, 5000)
Time_rec3 <- matrix(NA,nrow=length(r.ss),ncol = 4)
colnames(Time_rec3) <- c("mean", "median", "sd", "max")
difn_result <- list()
count1 <- 1
for(r in r.ss){
  cputime <- vector()
  theta_all <- vector()
  for (i in 1:nrep){
    cat(10000+r+i)
    cat(" ")
    #load("SUSY_pre.Rdata")
    #rm("testdata")
    #gc()
    set.seed(30231102+i)
    ## pilot sampling
    NDATA <- N
    index0 <- sample(NDATA, r0, replace = FALSE)
    # step 1 data
    data0 <- DATA[index0,]
    X0 <- data0[,-1]
    Y0 <- data0[, 1]
    # step 2 data
    data <- DATA[-index0,]
    dataN <- nrow(data)
    #rm("DATA")
    gc()
    # OSMAC_L
    set.seed(20231112+i)
    
    pilot_theta <- as.vector(getMLE(X0, Y0)$par)
    P0 <- 1-1/(1+exp(X0%*%pilot_theta))
    Nphi0 <- mean(abs(Y0-P0) * sqrt(rowSums(X0^2)))*dataN
    
    
    #Possion sampling
    init.time <- Sys.time()
    sample.OSMAC <- OSMAC(x = data, r = r, psi = Nphi0, theta = pilot_theta)
    end.time <- Sys.time()
    cputime <- c(cputime, end.time-init.time)
    #cputime
    #tail(sample.OSMAC,1)[2]
    
    X <- rbind(sample.OSMAC[,-c(1,p+3)], X0)
    Y <- c(sample.OSMAC[,1], Y0)
    weight <- c(1/sample.OSMAC[,p+3], rep(1,r0))/(dataN + r0)
    theta_new <- getMLE(X ,Y, weight)$par
    end.time <- Sys.time()
    theta_all <- cbind(theta_all, theta_new)
    #rm("data")
  }
  Time_rec3[count1,1:4] <- c(mean(cputime),median(cputime),sd(cputime),max(cputime)) 
  difn_result[[count1]] <- theta_all
  count1 <- count1 + 1
}
theta_true <- theta_logistic
k <- 1
for (j in 1:length(r.ss)) {
  bias <- (difn_result[[j]]-theta_true)^2
  
  MSE <- colSums(bias)
  EMSE[j,k] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  sdMSE[j,k] <- sd(MSE, na.rm = T)
}
###########################################

################ RBIMP ###################
#r.ss <- c(2000, 5000)
Time_rec <- matrix(NA,nrow=length(r.ss),ncol = 4)
colnames(Time_rec) <- c("mean", "median", "sd", "max")
difn_result <- list()
count1 <- 1
for(r in r.ss){
  cputime <- vector()
  theta_all <- vector()
  for (i in 1:nrep){
    set.seed(30231102+i)
    cat(20000+r+i)
    cat(" ")
    #load("SUSY_pre.Rdata")
    #rm("testdata")
    #gc()
    NDATA <- N
    index0 <- sample(NDATA, r0, replace = FALSE)
    # step 1 data
    data0 <- DATA[index0,]
    X0 <- data0[,-1]
    Y0 <- data0[, 1]
    data <- DATA[-index0,]
    dataN <- nrow(data)
    #rm("DATA")
    gc()
    # pilot
    pilot_theta <- as.vector(getMLE(X0, Y0)$par)
    P0 <- 1-1/(1+exp(X0%*%pilot_theta))
    Nphi0 <- mean(abs(Y0-P0) * sqrt(rowSums(X0^2)))*dataN
    
    
    #Possion sampling
    init.time <- Sys.time()
    sample.RBIMP <- RBIMP(x = data, r = r, psi = Nphi0, theta = pilot_theta)
    end.time <- Sys.time()
    cputime <- c(cputime, end.time-init.time)
    
    rr <- (nrow(sample.RBIMP)-3)/2
    X <- rbind(sample.RBIMP[1:(rr+2),-c(1,p+3)], X0)
    Y <- c(sample.RBIMP[1:(rr+2),1], Y0)
    weight <- 1/sample.RBIMP[1:rr,p+3]
    Z2 <- tail(sample.RBIMP,rr)
    Z1 <- sample.RBIMP[rr+3,]
    wZ2 <- Z2*weight
    
    difference_mean <- Z1-colSums(wZ2)
    
    Z2_seondmoment <- t(Z2)%*%(wZ2)
    
    weight_p <- pmax(0,c(1 + t(difference_mean)%*%solve(Z2_seondmoment, t(Z2))))
    weight <- weight*weight_p
    
    weight <- c(weight, sample.RBIMP[(rr+1):(rr+2),p+3], rep(1,r0))/(dataN + r0)
    theta_new <- getMLE(X ,Y, weight)$par
    
    theta_all <- cbind(theta_all, theta_new)
    #rm("data")
  }
  Time_rec[count1,1:4] <- c(mean(cputime),median(cputime),sd(cputime),max(cputime)) 
  difn_result[[count1]] <- theta_all
  count1 <- count1 + 1
}
k <- 2
for (j in 1:length(r.ss)) {
  bias <- (difn_result[[j]]-theta_true)^2
  
  MSE <- colSums(bias)
  EMSE[j,k] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  sdMSE[j,k] <- sd(MSE, na.rm = T)
  
}
#########################################

rm("DATA")
rm("data")
#######################################
save(list = ls(), file=paste0("SUSY_logistic_CompTime_p=",p,".RData"))
