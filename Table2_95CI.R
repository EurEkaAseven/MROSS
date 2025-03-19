rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(MASS)
library(matrixStats)
library(doParallel) 
library(foreach)
source("basic/DataGenerator.R")
source("basic/PA.R")
source("basic/getMLE.R")

packs <- c("mvtnorm","MASS","matrixStats")
CI_rep <- function(i){
  epsilon <- Cs
  set.seed(30231102+i)
    DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
  ## pilot sampling
  NDATA <- nrow(DATA)
  index0 <- sample(NDATA, r0, replace = FALSE)
  # step 1 data
  data0 <- DATA[index0,]
  X0 <- data0[,-1]
  Y0 <- data0[, 1]
  # step 2 data
  data <- DATA[-index0,]
  dataX <- data[,-1]
  dataY <- data[, 1]
  dataN <- nrow(data)
  pilot_theta <- as.vector(getMLE(X0, Y0)$par)
    ############ Uniform ##############
    set.seed(20231112+i)
      Pi.unif <- rep(r/dataN, dataN)
      
      sample_tmp <- (1:dataN)[runif(dataN) <= Pi.unif]
      X1 <- dataX[sample_tmp,]
      X <- rbind(X1, X0)
      Y <- c(dataY[sample_tmp], Y0)
      Pi.unif <- Pi.unif[sample_tmp]
      weight <- c(1/Pi.unif, rep(1,r0))
      fit.unif <- getMLE(X, Y, weight)
      theta.unif <- fit.unif$par
      
      P0 <- 1-1/(1+exp(X0%*%theta.unif))
      P1 <- 1-1/(1+exp(X1%*%theta.unif))
      H <- fit.unif$Infor
      Vc <- t(X0)%*%diag(c((Y0-P0)^2))%*%X0
      Vc <- (t(X1)%*%diag(c((dataY[sample_tmp]-P1)^2*(Pi.unif^(-2))))%*%X1-t(X1)%*%diag(c((dataY[sample_tmp]-P1)^2/Pi.unif))%*%X1)
      sigma.unif <- diag(solve(H)%*%Vc%*%solve(H))[ind.beta]

    ############# OSMAC #################
    set.seed(20231112+i)
      distant <- dataX%*%pilot_theta
      expcomp <- exp(distant)
      P <- 1-1/(1+expcomp)
      PmVc_tmp <- abs(dataY-P) * sqrt(rowSums(dataX^2))
      PmVc <- pmin(r*PmVc_tmp/sum(PmVc_tmp),1)
      
      #Possion sampling
      sample_tmp <- (1:dataN)[runif(dataN) <= PmVc]
      X1 <- dataX[sample_tmp,]
      X <- rbind(data[sample_tmp,-1], X0)
      Y <- c(data[sample_tmp, 1], Y0)
      PmVc <- PmVc[sample_tmp]
      weight <- c(1/PmVc, rep(1,r0))
      fit.mVc <- getMLE(X, Y, weight)
      theta.mVc <- fit.mVc$par
      
      P0 <- 1-1/(1+exp(X0%*%theta.mVc))
      P1 <- 1-1/(1+exp(X1%*%theta.mVc))
      H <- fit.mVc$Infor
      Vc <- t(X0)%*%diag(c((Y0-P0)^2))%*%X0
      Vc <- Vc + (t(X1)%*%diag(c((dataY[sample_tmp]-P1)^2*(PmVc^(-2))))%*%X1-t(X1)%*%diag(c((dataY[sample_tmp]-P1)^2/PmVc))%*%X1)
      
      sigma.mVc <- diag(solve(H)%*%Vc%*%solve(H))[ind.beta]
      
    ##############  RBIMP ############### 
    set.seed(20231112+i)
      distant <- dataX%*%pilot_theta
      expcomp <- exp(distant)
      P <- 1-1/(1+expcomp)
      ind_plus <- (1:dataN)[distant > epsilon & dataY == 1]
      ind_minu <- (1:dataN)[distant < (-epsilon) & dataY == 0]
      npos <- length(ind_plus)
      nneg <- length(ind_minu)
      if (npos>1){
        xbarpos <- colMeans(dataX[ind_plus,])
      }else{xbarpos <- vector()}
      if (nneg>1){
        xbarneg <- colMeans(dataX[ind_minu,])
      }else{xbarneg <- vector()}
      tdataN <- dataN-npos-nneg
      {if (nneg>0){
        if (npos>0){ind_s <- (1:dataN)[-c(ind_minu,ind_plus)]}
        else{ind_s <- (1:dataN)[-ind_minu]}
      }else{
        if (npos>0){ind_s <- (1:dataN)[-ind_plus]}
        else{ind_s <- (1:dataN)}
      }
      }
      if(r >= tdataN){
        X <- data[ind_s,-1]
        Y <- data[ind_s, 1]
        weight <- rep(1,length(ind_s))
      }else{
        ############################## Projection ##################################
        expcomp <- expcomp[ind_s]
        P <- 1-1/(1+expcomp)
        gradient <- (dataY[ind_s]-P)*dataX[ind_s,]
        
        Pi_tmp <- sqrt(rowSums(gradient^2))
        Pi <- pmin(1,r*Pi_tmp/sum(Pi_tmp))
        
        subdata_ind <- (runif(tdataN) <= Pi)
        sample_tmp <- (ind_s)[subdata_ind]
        X <- dataX[sample_tmp,]
        Y <- dataY[sample_tmp]
        weight <- 1/Pi[subdata_ind]
        
        Z1 <- c(tdataN, sum(dataY[ind_s]), colSums(gradient))
        Z2 <- cbind(1,Y,gradient[subdata_ind,])
        wZ2 <- Z2*weight
        
        difference_mean <- Z1-colSums(wZ2)
        Z2_seondmoment <- t(Z2)%*%wZ2
        weight_p <- pmax(0,as.vector(1 + t(difference_mean)%*%solve(Z2_seondmoment, t(Z2))))
        weight1 <- weight
        weight <- weight*weight_p
        X1 <- X
        Y1 <- Y
      }
      ############################## Combine #####################################
    
    if(npos > 1){
        X <- rbind(X,xbarpos)
        Y <- c(Y,1)
        weight <- c(weight, npos)
      }
      if(nneg > 1){
        X <- rbind(X,xbarneg)
        Y <- c(Y,0)
        weight <- c(weight, nneg)
      }
      
      X <- rbind(X, X0)
      Y <- c(Y, Y0)
      weight <- c(weight, rep(1,r0))
      
      fit.Imp <- getMLE(X ,Y, weight)
      theta.Imp <- fit.Imp$par
      
      H <- fit.mVc$Infor/(r0+dataN)
      P0 <- 1-1/(1+exp(X0%*%theta.Imp))
      P1 <- 1-1/(1+exp(X1%*%theta.Imp))
      score <- c(Y1-P1)*X1
      B <- solve(Z2_seondmoment, t(Z2)) %*% (weight1*score)
      aleph <- (score - Z2%*%B)
      Pi <- Pi[subdata_ind]

      Vc <- t(X0)%*%diag(c((Y0-P0)^2))%*%X0
      Vc <- Vc + (t(aleph)%*%diag((Pi^(-2)))%*%aleph) -t(aleph)%*%diag(1/(Pi))%*%aleph
      Vc <- Vc/length(r0 + ind_s)^2

      sigma.Imp <- diag(solve(H)%*%Vc%*%solve(H))[ind.beta]
  
      ############# MoreEff ###################
      r1 <- r0
      r2 <- r
      n <- dataN
      X <- dataX
      P  <- 1 - 1 / (1 + exp(X %*% pilot_theta))
      Y  <- dataY
      n1 <- sum(Y)
      n0 <- n - n1
      d <- ncol(X)
      
      ## First step
      PI.prop <- rep(1/(2*n0), n)
      PI.prop[Y==1] <- 1/(2*n1)
      idx.prop <- sample(1:n, r1, T, PI.prop)
      x.prop <- X[idx.prop,]
      y.prop <- Y[idx.prop]
      pinv.prop <- 1/PI.prop[idx.prop]
      fit.prop <- getMLE(x=x.prop, y=y.prop, w=1)
      beta.prop <- fit.prop$par - c(log(n0/n1), rep(0, d-1))
      P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
      
      if (is.na(beta.prop[1])) {
        warning("first stage not converge")
        next
      }
      
      p.propT <- 1 - 1 / (1 + exp(c(x.prop %*% fit.prop$par)))
      phi.propT <- p.propT * (1 - p.propT)
      ldd.prop <- t(x.prop) %*% (x.prop * phi.propT)
      
      psi.propT <- (y.prop - p.propT)^2
      psidd.prop <- t(x.prop) %*% (x.prop * psi.propT)
      
      p.prop <- P.prop[idx.prop]
      phi.prop <- p.prop * (1 - p.prop)
      ## ldd.prop <- t(x.prop) %*% (x.prop * phi.prop)
      
      ## mVc
      PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
      ## PI.mVc <- PI.mVc / sum(PI.mVc)
      idx.mVc <- sample(1:n, r2, T, PI.mVc)
      
      ## POImVcUW
      PI.poimVcuw <- r2 * PI.mVc / mean(PI.mVc[idx.prop] * pinv.prop)
      ## PI.poimVcuw <- r2 * PI.mVc
      idx.poimVcuw <- (1:n)[runif(n) <= PI.poimVcuw]
      x.poimVcuw <- X[c(idx.poimVcuw),]
      y.poimVcuw <- Y[c(idx.poimVcuw)]
      w.poimVcuw <- pmax(PI.poimVcuw[idx.poimVcuw], 1)
      fit.poimVcuw <- getMLE(x=x.poimVcuw, y=y.poimVcuw, w=w.poimVcuw)
      beta.poimVcuw <- fit.poimVcuw$par+beta.prop
      p.poimVcuw  <- 1 - 1 / (1 + exp(c(x.poimVcuw %*% fit.poimVcuw$par)))
      ## p.poimVcuw  <- 1 - 1 / (1 + exp(c(x.poimVcuw %*% beta.poimVcuw)))
      phi.poimVcuw <- p.poimVcuw * (1 - p.poimVcuw)
      ldd.poimVcuw <- t(x.poimVcuw) %*% (x.poimVcuw * phi.poimVcuw)
      beta.poimVcuwcb <- solve(ldd.prop+ldd.poimVcuw,
                               ldd.prop %*% beta.prop + ldd.poimVcuw %*% beta.poimVcuw)
      if (any(fit.poimVcuw$message != "Successful convergence")) {
        warning(paste("r2", r2, "not converge-.poimVcuw",
                      fit.poimVcuw$msg))
        beta.poimVcuwcb <- rep(NA, d)
      }
      ## V.poimVcuwcb <- solve(ldd.prop+ldd.poimVcuw)
      psi.poimVcuw <- (y.poimVcuw - p.poimVcuw)^2
      psidd.poimVcuw <- t(x.poimVcuw) %*% (x.poimVcuw * psi.poimVcuw)
      M.poimVcuwcb <- solve(ldd.prop+ldd.poimVcuw)
      V.poimVcuwcb <- M.poimVcuwcb %*% (psidd.prop+psidd.poimVcuw) %*% M.poimVcuwcb
      sigma.MEff <- diag(V.poimVcuwcb)[ind.beta]
      
  return(c(sigma.unif,sigma.mVc,sigma.MEff,sigma.Imp,
           theta.unif[ind.beta],theta.mVc[ind.beta],beta.poimVcuwcb[ind.beta],theta.Imp[ind.beta]))
}

# basic settings
ind.beta <- 2              # ind.beta = 1 for intercept term
nrep <-500                 # Replication size
N <- 5e5                   # Full data size
p <- 20                    # Dimension of X
r0 <- 1000                 # Step 1 sampling size
rs <- r0*seq(2,5,1)        # Step 2 sampling size

# Working model
wmodel <- "logistic"
# Generating model
gmodel <- "logistic" 
option <- "i" 

set.seed(20231102)
theta_true <- c(0.5,rep(0.5,p))      # True theta of logistic generating model

testdata <- GenerateData(N = 0.1*N, p = p, theta_true, 
                         option = option, model = gmodel)

if (wmodel == "logistic"){
  Cs <-3*log(10)
  compare_method <- c("Unif", "OSMAC_L", "MoreEff_Wang_L")
  if(gmodel != "logistic"){
    bigDATA <- GenerateData(N = 10*N, p = p, theta_true, option = option, model = gmodel)
    theta_true <- as.vector(getMLE(bigDATA[,-1], bigDATA[,1])$par)
  }
  PA_true <- PreAccuracy(testdata, theta_true, model = "logistic")
}

nm <- length(compare_method)

difn_result <- list()
for (j in 1:length(rs)){
  r <- rs[j]
  system.time({
    cl<- makeCluster(detectCores())  
    registerDoParallel(cl) 
    rep_result <- foreach(
      i=1:nrep,          
      .combine="rbind",  
      #.errorhandling = "pass",
      .packages=packs
    ) %dopar% {CI_rep(i)}
    stopCluster(cl)
  })
  difn_result[[j]] <- rep_result
}


ntheta <- nm + length(Cs) 
coverrate <- matrix(NA,nrow = length(rs),ncol = ntheta)
coverrate <- as.data.frame(coverrate)
rownames(coverrate) <- rs
colnames(coverrate)[1:nm] <- compare_method
colnames(coverrate)[(nm+1):(nm+length(Cs))] <- "Alg.1"
coverlen <- coverrate

for (j in 1:length(rs)) {
  sigma <- (difn_result[[j]])[,1:4]
  beta <- (difn_result[[j]])[,5:8]
  width <- 1.96*sqrt(sigma)
  coverlen[j,] <- colMeans(2*width)
  coverrate[j,] <- colMeans((theta_true[ind.beta]>=beta-width & theta_true[ind.beta]<=beta+width))
}
coverlen
coverrate

# save the results
if(wmodel == "logistic" & gmodel != "logistic"){rm("bigDATA")}
if(wmodel == "dwd" & gmodel != "Bdwd"){rm("bigDATA")}                   
save(list = ls(),
     file=paste0("beta_",ind.beta,"_95CI_wmodel=",as.character(wmodel),
                 "_gmodel=",as.character(gmodel),
                 "_option=",as.character(option),
                 "_p=", p, "_r0=", r0, "optimal=L",'.RData'))