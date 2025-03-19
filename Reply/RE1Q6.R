rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(MASS)
library(Rcpp)
source("basic/DataGenerator.R")
source("basic/getMLE.R")
source("Reply/GDgetMLE.R")
source("Reply/IRLSgetMLE.R")
sourceCpp("Rcpp/OSMAC.cpp")

nrep <-200               # Replication size
N <- 5e5                 # Full data size
p <- 20                  # Dimension of X
r0 <- 1000               # Step 1 sampling size
r.ss <- c(2000, 5000)    # Step 2 sampling size

#parameters of a decayed learning rate for GD
learning_rate0 <- 4 
alpha <- 0.6

# Working model
wmodel <- "logistic"
# Generating model
gmodel <- "logistic"

option <- "i"                       # Distribution setting of X
theta_true <- c(0.5,rep(0.5,p))      # True theta of logistic generating model
METHODS <- c("Full-IRLS","OSMAC-IRLS","OSMAC-GD","Full-GD")
nm <- length(METHODS)
EMSE <- matrix(NA,nrow = length(r.ss),ncol = nm)
EMSE <- as.data.frame(EMSE)
rownames(EMSE) <- r.ss
colnames(EMSE) <- METHODS
sdMSE <- EMSE

############## OSMAC ################
Time_before <- data.frame(Method = rep(METHODS, length(r.ss)),
                         r = rep(r.ss, each = nm),
                         Mean = 0, Median = 0, sd = 0, Max = 0)
Time_samp <- Time_before
Time_solve <- Time_before
Time_whole <- Time_before
difn_result <- list()
count1 <- 1
for(r in r.ss){
  cputime <- vector()
  cputime_samp <- vector()
  cputime_1 <- vector()
  cputime_2 <- vector()
  cputime_3 <- vector()
  cputime_4 <- vector()
  theta_all <- vector()
  for (i in 1:nrep){
    set.seed(30231102+i)
    DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
    #S4
    {if (count1 == 2){
     init.time_4 <- Sys.time()
     theta_new4 <- GDgetMLE(DATA[,-1], DATA[,1], w = 1/N)$par
     end.time_4 <- Sys.time()
     {if(is.na(theta_new4[1])){theta_new4 <- rep(NA,p+1)}
      else{theta_new4 <- theta_new4/norm(theta_new4[-1], "2")*norm(theta_true[-1], "2")}}
     #S1
     init.time_1 <- Sys.time()
     theta_new1 <- IRLSgetMLE(DATA[,-1], DATA[,1], w = 1/N)$par
     end.time_1 <- Sys.time()
     if(is.na(theta_new1[1])){theta_new1 <- rep(NA,p+1)}}
     else{
     init.time_4 <- Sys.time()
     theta_new4 <- rep(0,p+1)
     end.time_4 <- Sys.time()
      init.time_1 <- Sys.time()
      theta_new1 <- rep(0,p+1)
      end.time_1 <- Sys.time()}}
    
    ## pilot sampling
    NDATA <- nrow(DATA)
    index0 <- sample(NDATA, r0, replace = FALSE)
    # step 1 data
    data0 <- DATA[index0,]
    X0 <- data0[,-1]
    Y0 <- data0[, 1]
    # step 2 data
    data <- DATA[-index0,]
    dataN <- nrow(data)
    rm("DATA")
    gc()
    # OSMAC_L
    set.seed(20231112+i)
    init.time <- Sys.time()
    pilot_theta <- as.vector(getMLE(X0, Y0)$par)
    P0 <- 1-1/(1+exp(X0%*%pilot_theta))
    Nphi0 <- mean(abs(Y0-P0) * sqrt(rowSums(X0^2)))*dataN
    
    
    
    #Possion sampling
    init.time_samp <- Sys.time()
    sample.OSMAC <- OSMAC(x = data, r = r, psi = Nphi0, theta = pilot_theta)
    end.time_samp <- Sys.time()
    
    X <- rbind(sample.OSMAC[,-c(1,p+3)], X0)
    Y <- c(sample.OSMAC[,1], Y0)
    weight <- c(1/sample.OSMAC[,p+3], rep(1,r0))/(dataN + r0)
    end.time <- Sys.time()
    
    
    #S2
    init.time_2 <- Sys.time()
    theta_new2 <- IRLSgetMLE(X ,Y, weight)$par
    end.time_2 <- Sys.time()
    i = i + 1
    if(is.na(theta_new2[1])){theta_new2 <- rep(NA,p+1)
    cat(i)}
    #S3
    init.time_3 <- Sys.time()
    theta_new3 <- GDgetMLE(X ,Y, weight)$par
    end.time_3 <- Sys.time()
    {if(is.na(theta_new3[1])){theta_new3 <- rep(NA,p+1)}
      else{theta_new3 <- theta_new3/norm(theta_new3[-1], "2")*norm(theta_true[-1], "2")}}
    
    cputime <- c(cputime, end.time-init.time)
    cputime_samp <- c(cputime_samp, end.time_samp-init.time_samp)
    cputime_1 <- c(cputime_1, end.time_1-init.time_1)
    cputime_2 <- c(cputime_2, end.time_2-init.time_2)
    cputime_3 <- c(cputime_3, end.time_3-init.time_3)
    cputime_4 <- c(cputime_4, end.time_4-init.time_4)
    theta_all <- cbind(theta_all, c(theta_new1, theta_new2, theta_new3, theta_new4))
    rm("data")
  }
  Time_before[((count1-1)*nm+2):(count1*nm-1), 3:6] <- matrix(rep(c(mean(cputime),median(cputime),sd(cputime),max(cputime)),each = nm-2), ncol = 4) 
  Time_samp[((count1-1)*nm+2):(count1*nm-1), 3:6] <- matrix(rep(c(mean(cputime_samp),median(cputime_samp),sd(cputime_samp),max(cputime_samp)),each = nm-2), ncol = 4) 
  Time_solve[((count1-1)*nm+1), 3:6] <- c(mean(cputime_1),median(cputime_1),sd(cputime_1),max(cputime_1))
  Time_solve[((count1-1)*nm+2), 3:6] <- c(mean(cputime_2),median(cputime_2),sd(cputime_2),max(cputime_2))
  Time_solve[((count1-1)*nm+3), 3:6] <- c(mean(cputime_3),median(cputime_3),sd(cputime_3),max(cputime_3))
  Time_solve[((count1-1)*nm+4), 3:6] <- c(mean(cputime_4),median(cputime_4),sd(cputime_4),max(cputime_4))
  
  Time_whole[((count1-1)*nm+1), 3:6] <- c(mean(cputime_1),median(cputime_1),
                                          sd(cputime_1),max(cputime_1))
  Time_whole[((count1-1)*nm+2), 3:6] <- c(mean(cputime_2+cputime),median(cputime_2+cputime),
                                          sd(cputime_2+cputime),max(cputime_2+cputime))
  Time_whole[((count1-1)*nm+3), 3:6] <- c(mean(cputime_3+cputime),median(cputime_3+cputime),
                                          sd(cputime_3+cputime),max(cputime_3+cputime))
  Time_whole[((count1-1)*nm+4), 3:6] <- c(mean(cputime_4),median(cputime_4),
                                          sd(cputime_4),max(cputime_4))
  difn_result[[count1]] <- theta_all
  count1 <- count1 + 1
}

for (j in 1:length(r.ss)) {
  bias <- (difn_result[[j]]-theta_true)^2
  for (k in 1:nm) {
    MSE <- colSums(bias[seq((k-1)*(p+1)+1,k*(p+1),1),])
    EMSE[j,k] <- mean(MSE, na.rm = T)
    cat(sum(is.na(MSE)))
    cat(" ")
    sdMSE[j,k] <- sd(MSE, na.rm = T)
  }
}
###########################################


#######################################
save(list = ls(), file=paste0("SampleTime_p=",p,".RData"))
