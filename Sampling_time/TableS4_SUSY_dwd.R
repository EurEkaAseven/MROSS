rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Reply1")
#################################################
library(MASS)
library(Rcpp)
library(kerndwd)
source("basic/DataGenerator.R")
sourceCpp("Rcpp/OSMAC_dwd.cpp")
sourceCpp("Rcpp/RBIMP_dwd.cpp")

load("realdata/SUSY_pre.Rdata")
DATA[,1] <- 2*DATA[,1]-1
rm("testdata")
gc()

nrep <-200                  # Replication size
r0 <- 1000                   # Step 1 sampling size
r.ss <- c(2000,5000)
# Working model
wmodel <- "dwd"
lambda = 0
kern = vanilladot()

theta_full <- theta_dwd

EMSE <- matrix(NA,nrow = length(r.ss),ncol = 2)
EMSE <- as.data.frame(EMSE)
rownames(EMSE) <- r.ss
colnames(EMSE) <- c("OSMAC","RBIMP")
sdMSE <- EMSE

############## OSMAC ###################
#r.ss <- c(2000,5000)
Time_OSMAC <- matrix(NA,nrow=length(r.ss),ncol = 4)
Time_OSMAC_samp <- Time_OSMAC
difn_result <- list()
count1 <- 1
for(r in r.ss){
  cputime <- vector()
  cputime_samp <- vector()
  theta_all <- vector()
  for (i in 1:nrep){
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
    init.time <- Sys.time()
    pilot_theta <- as.vector(kerndwd(X0[,-1], Y0, 
                                     kern = kern,
                                     qval = 1, lambda = lambda,
                                     eps = 1e-5, maxit = 1e+5)$alpha)
    u0 <- (X0%*%pilot_theta)*Y0
    Nphi0 <- mean(ifelse(u0>0.5, 1/u0^2/4, 1) * sqrt(rowSums(X0^2)))*dataN
    
    
    #Possion sampling
    init.time_samp <- Sys.time()
    sample.OSMAC <- OSMAC_dwd(x = data, r = r, psi = Nphi0, theta = pilot_theta)
    end.time_samp <- Sys.time()
    
    X <- rbind(sample.OSMAC[,-c(1,p+3)], X0)
    Y <- c(sample.OSMAC[,1], Y0)
    weight <- c(1/sample.OSMAC[,p+3], rep(1,r0))/(dataN + r0)
    theta_new <- kerndwd(X[,-1], Y, wt = weight,
                         kern = kern,
                         qval = 1, lambda = lambda,
                         eps = 1e-5, maxit = 1e+5)$alpha
    end.time <- Sys.time()
    cputime <- c(cputime, end.time-init.time)
    cputime_samp <- c(cputime_samp, end.time_samp-init.time_samp)
    theta_all <- cbind(theta_all, theta_new)
    rm("data")
  }
  Time_OSMAC[count1,1:4] <- c(mean(cputime),median(cputime),sd(cputime),max(cputime))
  Time_OSMAC_samp[count1,1:4] <- c(mean(cputime_samp),median(cputime_samp),
                                   sd(cputime_samp),max(cputime_samp)) 
  difn_result[[count1]] <- theta_all
  count1 <- count1 + 1
}
k <- 1
for (j in 1:length(r.ss)) {
  bias <- (difn_result[[j]]-theta_full)^2
  
  MSE <- colSums(bias)
  EMSE[j,k] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  sdMSE[j,k] <- sd(MSE, na.rm = T)
  
} 
###########################################

################ RBIMP ###################
#r.ss <- c(2000,5000)
Time_RBIMP <- matrix(NA,nrow=length(r.ss),ncol = 4)
Time_RBIMP_samp <- Time_OSMAC
difn_result <- list()
count1 <- 1
for(r in r.ss){
  cputime <- vector()
  cputime_samp <- vector()
  theta_all <- vector()
  for (i in 1:nrep){
    set.seed(30231102+i)
    #DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
    #DATA[,1] <- 2*DATA[,1]-1
    ## pilot sampling
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
    init.time <- Sys.time()
    pilot_theta <- as.vector(kerndwd(X0[,-1], Y0, 
                                     kern = kern,
                                     qval = 1, lambda = lambda,
                                     eps = 1e-5, maxit = 1e+5)$alpha)
    u0 <- (X0%*%pilot_theta)*Y0
    Nphi0 <- mean(ifelse(u0>0.5, 1/u0^2/4, 1) * sqrt(rowSums(X0^2)))*dataN
    
    
    #Possion sampling
    init.time_samp <- Sys.time()
    sample.RBIMP <- RBIMP_dwd(x = data, r = r, psi = Nphi0, theta = pilot_theta)
    end.time_samp <- Sys.time()
    
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
    na.index <- which(is.na(Y))
    if(length(na.index) > 0){
    X <- X[-na.index,]
    Y <- Y[-na.index]
    weight <- weight[-na.index]}
    
    theta_new <- kerndwd(X[,-1], Y, wt = weight,
                         kern = kern,
                         qval = 1, lambda = lambda,
                         eps = 1e-5, maxit = 1e+5)$alpha
    end.time <- Sys.time()
    cputime <- c(cputime, end.time-init.time)
    cputime_samp <- c(cputime_samp, end.time_samp-init.time_samp)
    theta_all <- cbind(theta_all, theta_new)
    rm("data")
  }
  Time_RBIMP[count1,1:4] <- c(mean(cputime),median(cputime),sd(cputime),max(cputime))
  Time_RBIMP_samp[count1,1:4] <- c(mean(cputime_samp),median(cputime_samp),
                                   sd(cputime_samp),max(cputime_samp)) 
  difn_result[[count1]] <- theta_all
  count1 <- count1 + 1
}
k <- 2
for (j in 1:length(r.ss)) {
  bias <- (difn_result[[j]]-theta_full)^2
  
  MSE <- colSums(bias)
  EMSE[j,k] <- mean(MSE, na.rm = T)
  cat(sum(is.na(MSE)))
  cat(" ")
  sdMSE[j,k] <- sd(MSE, na.rm = T)
  
}
#########################################

#######################################
rm("DATA")
rm("data")
save(list = ls(), file=paste0("dwdCompTime_p=",p,".RData"))
