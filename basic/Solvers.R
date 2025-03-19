solver <- function(data, data0, n, method, wmodel,
                   epsilon = 2, lambda = 0, validata = NULL,  delta = NULL){
  dataX <- data[,-1]
  dataY <- data[, 1]
  n0 <- nrow(data0)
  dataN <- nrow(data)
  
  if(wmodel == "logistic"){
    {
      {pilot_theta <- as.vector(getMLE(data0[,-1], data0[,1])$par)}
      
      if (method == "Unif"){
        Pi <- rep(n/dataN, dataN)
        
        sample_tmp <- sample_tmp <- (1:dataN)[runif(dataN) <= Pi]
        X <- rbind(data[sample_tmp,-1],data0[,-1])
        Y <- c(data[sample_tmp, 1],data0[,1])
        
        theta_new <- getMLE(X ,Y)$par
        result <- list(theta = as.vector(theta_new))
      }
      else if (method == "OSMAC_L"){
        ##################### Calculate PmVc ###########################
        expcomp <- exp(dataX%*%pilot_theta)
        P <- 1-1/(1+expcomp)
        PmVc_tmp <- abs(dataY-P) * sqrt(rowSums(dataX^2))
        PmVc <- pmin(n*PmVc_tmp/sum(PmVc_tmp),1)
        
        #Possion sampling
        sample_tmp <- (1:dataN)[runif(dataN) <= PmVc]
        X <- rbind(data[sample_tmp,-1], data0[,-1])
        Y <- c(data[sample_tmp, 1], data0[, 1])
        #weight <- c(n/(n+n0)/dataN/PmVc[sample_tmp], 1/(n+n0)/dataN/PI0)
        weight <- c(1/PmVc[sample_tmp], rep(1,n0))/(dataN + n0)
        
        theta_new <- getMLE(X ,Y, weight)$par
        result <- list(theta = as.vector(theta_new))
      }
      else if(method == "OSMAC_A"){
        ##################### Calculate PmVc ###########################
        P0 <-1-1/(1+exp(data0[,-1]%*%pilot_theta))
        P <- 1-1/(1+exp(dataX%*%pilot_theta))
        H <- t(data0[,-1])%*%(as.vector(P0 * (1-P0)) * data0[,-1])/nrow(data0)
        PmV_tmp <- abs(dataY-P) * sqrt(colSums(solve(H, t(dataX))^2))
        
        PmV <- pmin(n*PmV_tmp/sum(PmV_tmp),1)
        
        #Possion sampling
        sample_tmp <- (1:dataN)[runif(dataN) <= PmV]
        X <- rbind(data[sample_tmp,-1], data0[,-1])
        Y <- c(data[sample_tmp, 1], data0[, 1])
        #weight <- c(n/(n+n0)/dataN/PmV[sample_tmp], 1/(n+n0)/dataN/PI0)
          weight <- c(1/PmV[sample_tmp], rep(1,n0))/(dataN + n0)
        
        theta_new <- getMLE(X ,Y, weight)$par
        result <- list(theta = as.vector(theta_new))
      }
      else if(method == "MoreEff_Wang_A"){
        theta_new <- getMoreEff(rbind(dataX,data0[,-1]), c(dataY,data0[,1]),
                                n0, n, type = "PomMSE")$par
        result <- list(theta = as.vector(theta_new))
      }
      else if(method == "MoreEff_Wang_L"){
        theta_new <- getMoreEff(rbind(dataX,data0[,-1]), c(dataY,data0[,1]),
                                n0, n, type = "PomVc")$par
        result <- list(theta = as.vector(theta_new))
      }
      else if(method == "IBOSS"){
        sdata <- samp.data(data = list(y = dataY, x.matrix = dataX),
                           beta = pilot_theta,
                           num.covariate = ncol(data)-1,
                           num.datalines = dataN,
                           num.samp = n, delta = delta)
        theta_IBOSS <- getMLE(sdata$x.matrix, sdata$y)$par
        result <- list(theta = as.vector(theta_IBOSS))
        
      }
      else if(method == "Improve_L"){
        ############################### Truncation #################################

        distant <- dataX%*%pilot_theta
        #index <- rep(0, dataN)
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
        #tdata <- data[-c(ind_minu,ind_plus),]
        tdataN <- dataN-npos-nneg
        {if (nneg>0){
                    if (npos>0){ind_s <- (1:dataN)[-c(ind_minu,ind_plus)]}
                    else{ind_s <- (1:dataN)[-ind_minu]}
                }else{
                    if (npos>0){ind_s <- (1:dataN)[-ind_plus]}
                    else{ind_s <- (1:dataN)}
                }
        }
        if(n >= tdataN){
          X <- data[ind_s,-1]
          Y <- data[ind_s, 1]
          weight <- rep(1,length(ind_s))
        }else{
          ############################## Projection ##################################
          # tic()
          # tX <- tdata[,-1]
          # tY <- tdata[, 1]
          # cat("index")
          # toc()
            
          expcomp <- exp(distant[ind_s])
          P <- 1-1/(1+expcomp)
          #phi <- (Y[ind_s]-P)

          gradient <- (dataY[ind_s]-P)*dataX[ind_s,]

          Pi_tmp <- sqrt(rowSums(gradient^2))
          #Pi_tmp <- phi * sqrt(rowSums(tX^2))
          Pi <- pmin(1,n*Pi_tmp/sum(Pi_tmp))

          
          #sample_tmp <- rbinom(tdataN,1,Pi)
          subdata_ind <- (runif(tdataN) <= Pi)
          sample_tmp <- (ind_s)[subdata_ind]
          X <- dataX[sample_tmp,]
          Y <- dataY[sample_tmp]
          weight <- 1/Pi[subdata_ind]

          Z1 <- c(tdataN, sum(dataY[ind_s]), colSums(gradient))

          Z2 <- cbind(1,Y,gradient[subdata_ind,])
          wZ2 <- Z2*weight
          
          difference_mean <- (Z1)-colSums(wZ2)
          Z2_seondmoment <- t(Z2)%*%(wZ2)
          weight_p <- pmax(0,as.vector(1 + t(difference_mean)%*%solve(Z2_seondmoment, t(Z2))))
          weight <- weight*weight_p
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
        
        X <- rbind(X, data0[,-1])
        Y <- c(Y, data0[, 1])
        #weight <- c(n/(n+n0)/dataN*weight, rep(1/(n+n0), n0))
        #weight <- c(n/(n+n0)/dataN*weight, 1/(n+n0)/dataN/PI0)
          weight <- c(weight, rep(1,n0))/(dataN + n0)
        
        theta_new <- getMLE(X ,Y, weight)$par

        result <- list(theta = as.vector(theta_new))}
      else if(method == "Improve_A"){#Mainly for plotting figure 2 here 
        ############################### Truncation #################################
        distant <- dataX%*%pilot_theta
        index <- rep(0, dataN)
        for (i in 1:dataN) {
          if(distant[i] > epsilon & dataY[i] == 1){ index[i] <- 2 }
          else if(distant[i] < -epsilon & dataY[i] == 0){ index[i] <- 3 }
        }
        npos <- sum(index==2)
        nneg <- sum(index==3)
        cat(npos)
        cat(" ")
        if (npos>1){
          xbarpos <- colMeans(data[index==2,-1])
        }else{xbarpos <- vector()}
        if (nneg>1){
          xbarneg <- colMeans(data[index==3,-1])
        }else{xbarneg <- vector()}
        tdata <- data[index==0,]
        
        if(n >= nrow(tdata)){
          X <- tdata[,-1]
          Y <- tdata[, 1]
          weight <- rep(1,nrow(tdata))
          
        }else{
          ############################## Projection ##################################
          tX <- tdata[,-1]
          tY <- tdata[, 1]
          
          #P0 <-1-1/(1+exp(data0[,-1]%*%pilot_theta))
          P <- 1-1/(1+exp(tX%*%pilot_theta))
          H <- t(tX)%*%(as.vector(P * (1-P)) * tX) / nrow(tX)
          #H <- t(data0[,-1])%*%(as.vector(P0 * (1-P0)) * data0[,-1]) / nrow(data0)
          PmV_tmp <- abs(tY-P) * sqrt(colSums(solve(H, t(tX))^2))
          PmV <- pmin(n*PmV_tmp/sum(PmV_tmp),1)
          
          Projection <- getProjection(tdata, n = n, 
                                      pilot_theta = pilot_theta, 
                                      Pi = PmV)
          weight <- Projection$weight
          weight_p <- Projection$weight_p
          X <- Projection$sX
          Y <- Projection$sY
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
        
        X <- rbind(X, data0[,-1])
        Y <- c(Y, data0[, 1])
        weight <- c(n/(n+n0)/dataN*weight, rep(1/(n+n0), n0))
        
        theta_new <- getMLE(X ,Y, weight)$par
        result <- list(np = npos, nn = nneg, xbar = cbind(xbarpos,xbarneg),
                       theta = as.vector(theta_new), 
                       weight = weight[1:length(weight_p)], weight_p = weight_p, 
                       sdata = cbind(Y,X)[1:length(weight_p),],
                       tdata = tdata)
      }
      else{cat("Invalid method.")}
    }
  }
  else if(wmodel == "dwd"){
    pilot_theta <- as.vector(getDWD(data0[,-1], data0[,1],
                                    lambda = lambda, validata = validata)$par)
    
    if(method == "Lopt"){
      u <- (dataX%*%pilot_theta)*(2*dataY-1)
      derivation <- ifelse(u<=0.5, 1, 1/u^2/4)
      Pi_tmp <- derivation * sqrt(rowSums(dataX^2))
      Pi <- pmin(1,n*Pi_tmp/sum(Pi_tmp))
      
      sample_tmp <- (1:dataN)[runif(dataN) <= Pi]
      X <- rbind(data[sample_tmp,-1], data0[,-1])
      Y <- c(data[sample_tmp, 1], data0[, 1])
      #weight <- c(n/(n+n0)/dataN/Pi[sample_tmp], 1/(n+n0)/dataN/PI0)
      weight <- c(1/Pi[sample_tmp], rep(1,n0))/(dataN + n0)
        
      theta_new <- getDWD(X, Y, weight = weight, lambda = lambda, validata = validata)$par
      result <- list(theta = as.vector(theta_new))
    }
    else if(method == "Aopt"){
      #pilot_theta <- theta_true
        
      u <- (dataX%*%pilot_theta)*(2*dataY-1)
      u0 <- data0[,-1]%*%pilot_theta*(2*data0[,1]-1)
      derivation <- ifelse(u<=0.5, 1, 1/u^2/4)
      sec_derivation0 <- ifelse(u0<=0.5, 0, 1/u0^3/2)
      H <- t(data0[,-1])%*%(as.vector(sec_derivation0) * data0[,-1]) / nrow(data0)
      #sec_derivation <- ifelse(u<=0.5, 0, 1/u^3/2)
      #H <- t(dataX)%*%(as.vector(sec_derivation) * dataX) / dataN
        
      PmV_tmp <- derivation * sqrt(colSums(solve(H, t(dataX))^2))
      PmV <- pmin(1,n*PmV_tmp/sum(PmV_tmp))
      
      #shrink_rate <- 0.1
      #PmV <- (1-shrink_rate)*PmV + shrink_rate*n/dataN
      
      sample_tmp <- (1:dataN)[runif(dataN) <= PmV]
      X <- rbind(data[sample_tmp,-1], data0[,-1])
      Y <- c(data[sample_tmp, 1], data0[, 1])
      #weight <- c(n/(n+n0)/dataN/PmV[sample_tmp], 1/(n+n0)/dataN/PI0)
      weight <- c(1/PmV[sample_tmp], rep(1,n0))/(dataN + n0)
      
      theta_new <- getDWD(X, Y, weight = weight, lambda = lambda, validata = validata)$par
      result <- list(theta = as.vector(theta_new))
    }
    else if(method == "Unif"){
        Pi <- rep(n/dataN, dataN)
        
        sample_tmp <- sample_tmp <- (1:dataN)[runif(dataN) <= Pi]
        X <- rbind(data[sample_tmp,-1],data0[,-1])
        Y <- c(data[sample_tmp, 1],data0[,1])
      
      theta_new <- getDWD(X, Y, lambda = lambda, validata = validata)$par
      result <- list(theta = as.vector(theta_new))
    }
    else if(method == "Improve_L"){
      ############################### Truncation #################################
      distant <- dataX%*%pilot_theta
      #index <- rep(0, dataN)
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
      #tdata <- data[-c(ind_minu,ind_plus),]
      tdataN <- dataN-npos-nneg
      {if (nneg>0){
        if (npos>0){ind_s <- (1:dataN)[-c(ind_minu,ind_plus)]}
        else{ind_s <- (1:dataN)[-ind_minu]}
      }else{
        if (npos>0){ind_s <- (1:dataN)[-ind_plus]}
        else{ind_s <- (1:dataN)}
      }
      }
      if(n >= tdataN){
        X <- data[ind_s,-1]
        Y <- data[ind_s, 1]
        weight <- rep(1,length(ind_s))
      }else{
        ############################## Projection ##################################
        u <- distant[ind_s]*(2*dataY[ind_s]-1)
        gradient <- c((2*dataY[ind_s]-1)*ifelse(u>0.5, -1/u^2/4, -1))*dataX[ind_s,]
        
        Pi_tmp <- sqrt(rowSums(gradient^2))
        #Pi_tmp <- phi * sqrt(rowSums(tX^2))
        Pi <- pmin(1,n*Pi_tmp/sum(Pi_tmp))
        
        
        #sample_tmp <- rbinom(tdataN,1,Pi)
        subdata_ind <- (runif(tdataN) <= Pi)
        sample_tmp <- (ind_s)[subdata_ind]
        X <- dataX[sample_tmp,]
        Y <- dataY[sample_tmp]
        weight <- 1/Pi[subdata_ind]
        
        Z1 <- c(tdataN, sum(dataY[ind_s]), colSums(gradient))
        
        Z2 <- cbind(1,Y,gradient[subdata_ind,])
        wZ2 <- Z2*weight
        
        difference_mean <- (Z1)-colSums(wZ2)
        Z2_seondmoment <- t(Z2)%*%(weight*Z2)
        weight_p <- pmax(0,as.vector(1 + t(difference_mean)%*%solve(Z2_seondmoment, t(Z2))))
        weight <- weight*weight_p
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
      
      X <- rbind(X, data0[,-1])
      Y <- c(Y, data0[, 1])
      #weight <- c(n/(n+n0)/dataN*weight, rep(1/(n+n0), n0))
      #weight <- c(n/(n+n0)/dataN*weight, 1/(n+n0)/dataN/PI0)
      weight <- c(weight, rep(1,n0))/(dataN + n0)
      
      theta_new <- getDWD(X, Y, weight, lambda, validata)$par
      result <- list(theta = as.vector(theta_new))}
    else{cat("Invalid method.")}
    }
  else{cat("Invalid working model.")}
  
  return(result)
}

getDWD <- function(X, Y, weight = NULL, lambda, validata = NULL){
  if(is.null(weight)){
    DWD <- kerndwd(X[,-1], Y, 
                   kern = vanilladot(),
                   qval = 1, lambda = lambda,
                   eps = 1e-5, maxit = 1e+5)
  }else{
    DWD <- kerndwd(X[,-1], Y, 
                   kern = vanilladot(),
                   qval = 1, lambda = lambda,
                   wt = weight,
                   eps = 1e-5, maxit = 1e+5)
  }
  beta <- as.vector(DWD$alpha)
  
  if(is.null(validata)){
    result <- list(par = beta, lambda = lambda)
  }else{
    vY <- validata[, 1]
    vX <- validata[,-1]
    err <- nrow(X)
    for (k in 1:length(lambda)){
      pred <- vY*((vX %*% beta[,k]))
      temp_err <- sum(pred < 0)
      #cat(c(temp_err/nrow(vX)," "))
      if (temp_err < err){
        err <- temp_err
        index_best_lambda <- k
      }
    }
    result <- list(par = beta[,index_best_lambda], 
                   lambda = lambda[index_best_lambda])
  }
  return(result)
}

getMoreEff <- function(X, Y, r1, r2, type = "PomVc"){

  d <- ncol(X)
  n <- nrow(X)
  n1 <- sum(Y)
  n0 <- n - n1
  
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
    
    p.propT <- 1 - 1 / (1 + exp(c(x.prop %*% fit.prop$par)))
    phi.propT <- p.propT * (1 - p.propT)
    ldd.prop <- t(x.prop) %*% (x.prop * phi.propT)
    
    psi.propT <- (y.prop - p.propT)^2
    psidd.prop <- t(x.prop) %*% (x.prop * psi.propT)
    
    p.prop <- P.prop[idx.prop]
    phi.prop <- p.prop * (1 - p.prop)
    ## ldd.prop <- t(x.prop) %*% (x.prop * phi.prop)
    
    if (is.na(beta.prop[1])) {
      warning("first stage not converge")
      next
    }
  
  if(type == "mMSE"){
    ## mMSE
    W.prop <- solve(t(x.prop) %*% (x.prop * as.vector(phi.prop * pinv.prop)))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    ## PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
    ## mMSEUW
    x.mMSEuw <- X[c(idx.mMSE),]
    y.mMSEuw <- Y[c(idx.mMSE)]
    fit.mMSEuw <- getMLE(x=x.mMSEuw, y=y.mMSEuw, w=1)
    beta.mMSEuw <- fit.mMSEuw$par+beta.prop
    p.mMSEuw  <- 1 - 1 / (1 + exp(c(x.mMSEuw %*% fit.mMSEuw$par)))
    phi.mMSEuw <- p.mMSEuw * (1 - p.mMSEuw)
    ldd.mMSEuw <- t(x.mMSEuw) %*% (x.mMSEuw * phi.mMSEuw)
    beta.mMSEuwcb <- solve(ldd.prop+ldd.mMSEuw,
                           ldd.prop %*% beta.prop + ldd.mMSEuw %*% beta.mMSEuw)
    if (any(fit.mMSEuw$message != "Successful convergence"))
      warning(paste("itr", rr, "r2", r2, "not converge-.mMSEuw", fit.mMSEuw$msg))
    result <- list(par = as.vector(beta.mMSEuwcb))
  }else if(type == "mVc"){
    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    ## PI.mVc <- PI.mVc / sum(PI.mVc)
    idx.mVc <- sample(1:n, r2, T, PI.mVc)
    ## mVcUW
    x.mVcuw <- X[c(idx.mVc),]
    y.mVcuw <- Y[c(idx.mVc)]
    fit.mVcuw <- getMLE(x=x.mVcuw, y=y.mVcuw, w=1)
    beta.mVcuw <- fit.mVcuw$par+beta.prop
    p.mVcuw  <- 1 - 1 / (1 + exp(c(x.mVcuw %*% fit.mVcuw$par)))
    ## p.mVcuw  <- 1 - 1 / (1 + exp(c(x.mVcuw %*% beta.mVcuw)))
    phi.mVcuw <- p.mVcuw * (1 - p.mVcuw)
    ldd.mVcuw <- t(x.mVcuw) %*% (x.mVcuw * phi.mVcuw)
    beta.mVcuwcb <- solve(ldd.prop+ldd.mVcuw,
                          ldd.prop %*% beta.prop + ldd.mVcuw %*% beta.mVcuw)
    if (any(fit.mVcuw$message != "Successful convergence"))
      warning(paste("itr", rr, "r2", r2, "not converge-.mVcuw", fit.mVcuw$msg))
    result <- list(par = as.vector(beta.mVcuwcb))
  }else if(type == "PomMSE"){
    ## mMSE
    W.prop <- solve(t(x.prop) %*% (x.prop * as.vector(phi.prop * pinv.prop)))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    ## POImMSEUW
    PI.poimMSEuw <- r2 * PI.mMSE / mean(PI.mMSE[idx.prop] * pinv.prop)
    ## PI.poimMSEuw <- r2 * PI.mMSE
    idx.poimMSEuw <- (1:n)[runif(n) <= PI.poimMSEuw]
    x.poimMSEuw <- X[c(idx.poimMSEuw),]
    y.poimMSEuw <- Y[c(idx.poimMSEuw)]
    w.poimMSEuw <- pmax(PI.poimMSEuw[idx.poimMSEuw], 1)
    fit.poimMSEuw <- getMLE(x=x.poimMSEuw, y=y.poimMSEuw, w=w.poimMSEuw)
    beta.poimMSEuw <- fit.poimMSEuw$par+beta.prop
    p.poimMSEuw  <- 1 - 1 / (1 + exp(c(x.poimMSEuw %*% fit.poimMSEuw$par)))
    ## p.poimMSEuw  <- 1 - 1 / (1 + exp(c(x.poimMSEuw %*% beta.poimMSEuw)))
    phi.poimMSEuw <- p.poimMSEuw * (1 - p.poimMSEuw)
    ldd.poimMSEuw <- t(x.poimMSEuw) %*% (x.poimMSEuw * phi.poimMSEuw)
    beta.poimMSEuwcb <- solve(ldd.prop+ldd.poimMSEuw,
                              ldd.prop %*% beta.prop + ldd.poimMSEuw %*% beta.poimMSEuw)
    if (any(fit.poimMSEuw$message != "Successful convergence")) {
      warning(paste("itr", rr, "r2", r2, "not converge-.poimMSEuw",
                    fit.poimMSEuw$msg))
      beta.poimMSEuwcb <- rep(NA, d)
    }
    result <- list(par = as.vector(beta.poimMSEuwcb))
  }else if(type == "PomVc"){
    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
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
      warning(paste("itr", rr, "r2", r2, "not converge-.poimVcuw",
                    fit.poimVcuw$msg))
      beta.poimVcuwcb <- rep(NA, d)
    }

    result <- list(par = as.vector(beta.poimVcuwcb))
  }
  
  return(result)
}

getProjection <- function(data, n, Pi, pilot_theta){
  dataN <- nrow(data)
  sample_tmp <- rbinom(dataN,1,Pi)
  X <- data[sample_tmp==1,-1]
  Y <- data[sample_tmp==1, 1]
  weight <- 1/Pi[sample_tmp==1]

    gradient <- as.vector((2*data[,1]-1)*abs(data[,1]-(1-1/(1+exp(data[,-1]%*%pilot_theta)))))*data[,-c(1)]
    Z1 <- cbind(1,data[,1],gradient)
    Z2 <- Z1[sample_tmp==1,]
    wZ2 <- Z2*weight
    
    difference_mean <- colSums(Z1)-colSums(wZ2)
    Z2_seondmoment <- t(Z2)%*%(wZ2)
    weight_p <- as.vector(1 + t(difference_mean)%*%solve(Z2_seondmoment, t(Z2)))
    weight_p[weight_p<0] <- 0
    weight <- weight*weight_p
  
  result <- list(alpha = alpha, weight = weight, weight_p = weight_p, 
                 sX = X, sY = Y)
  return(result)
}
    
min.fun <- function(x,num.covariate){
  value=1/(x^2*(exp(x)/((1+exp(x))^2))^(num.covariate+1))
  return(value)  
}

