GenerateData<- function(N,p,theta,option = "i",model){

######## Generate Y by different Models ########################################  
  if (model == "logistic"){
    Sigma <- matrix(rep(0.5, p^2), nrow = p, ncol = p) + diag(rep(0.5,p))
    if (option == "i"){# Case 1
        Sigma <- matrix(nrow = p, ncol = p)
        for(i in 1:p){
          for(j in 1:p){
            Sigma[i,j] <- 0.5^abs(i-j)
         }
        }
      X <- mvrnorm(N, rep(0,p), Sigma)
    }else if (option == "ii"){# Case 2
        Sigma2 <- matrix(nrow = p, ncol = p)
        for(i in 1:p){
          for(j in 1:p){
            Sigma2[i,j] <- 0.5^abs(i-j)
         }
        }
      P <- rbinom(N, 1, 0.5)
      X <- (P*mvrnorm(N, rep(0,p), Sigma) + (1-P)*mvrnorm(N, rep(0,p), Sigma2))
    }else if (option == "iii"){# Case S1
        Sigma <- matrix(nrow = p, ncol = p)
        for(i in 1:p){
          for(j in 1:p){
            Sigma[i,j] <- 0.5^abs(i-j)
         }
        }
      X <- rmvt(N, delta = rep(0,p), sigma = Sigma, df = 3)
    }else if (option == "vii"){# Fig 1, 2ab
      X <- mvrnorm(N, rep(0,p), sqrt(2)*Sigma)
    }else if (option == "viii"){# Fig 2cd
      X <- vector()
      for (i in 1:p) {
        X <- cbind(X, runif(N, -4, 4))
      }
    }else if (option == "ix"){# Fig 2cd
      X <- vector()
      for (i in 1:p) {
        X <- cbind(X, rbeta(N, 4, 5))
      }
      X <- 10*X-5
    }
    X <- cbind(1,X)
    expcomp <- exp(X%*%theta)
    P <- expcomp/(1+expcomp)
    P[is.na(P)] <- 1
    Y <- rbinom(N, 1, P)
  }
  else if (model == "X_given_Y"){
    Y <- c(rep(1,N/2), rep(0,N/2))
    if (option == "iv"){# Case S2
      Yrate <- 0.8
      NP <- N*Yrate
      NN <- N - NP
      Y <- c(rep(1,NP), rep(0,NN))
      u1 <- c(rep(0,p/2),rep(1,p/2))
      u2 <- c(rep(0,p/2),rep(-1,p/2))
      Sigma <- matrix(nrow = p, ncol = p)
      for(i in 1:p){
        for(j in 1:p){
          Sigma[i,j] <- 0.5^abs(i-j)
        }
      }
      X1 <- rmvt(NP, delta = u1, sigma = Sigma, df = 3)
      X2 <- rmvt(NN, delta = u2, sigma = Sigma, df = 3)
      X <- rbind(X1,X2)
    }else if (option == "v"){#Case 3
      Sigma <- matrix(nrow = p, ncol = p)
      for(i in 1:p){
        for(j in 1:p){
          Sigma[i,j] <- 0.5^abs(i-j)
        }
      }
      Y <- c(rep(1,N/2), rep(0,N/2))
      u11 <- c(rep(0,ceiling(p/2)),rep(1,floor(p/2)))
      u12 <- c(rep(-1,ceiling(p/2)),rep(2,floor(p/2)))
      u13 <- rep(-1, p)
      u21 <- c(rep(0,ceiling(p/2)),rep(-1,floor(p/2)))
      u22 <- c(rep(1,ceiling(p/2)),rep(-2,floor(p/2)))
      u23 <- c(rep(1,ceiling(p/2)),rep(2,floor(p/2)))
      Sigma <- matrix(nrow = p, ncol = p)
      for(i in 1:p){
        for(j in 1:p){
          Sigma[i,j] <- 0.5^abs(i-j)
        }
      }
      
      index <- sample(1:4, N/2, replace = TRUE)
      Xp <- mvrnorm(N/2, u11, Sigma)
      Xp[index == 3] <-  mvrnorm(N/2, u12, Sigma)[index == 3]
      Xp[index == 4] <-  mvrnorm(N/2, u13, Sigma)[index == 4]
      index <- sample(1:4, N/2, replace = TRUE)
      Xn <- mvrnorm(N/2, u21, Sigma)
      Xn[index == 3] <-  mvrnorm(N/2, u22, Sigma)[index == 3]
      Xn[index == 4] <-  mvrnorm(N/2, u23, Sigma)[index == 4]
      
      X <- rbind(Xp,Xn)
    }else if (option == "vi"){# Case 4
      Y <- c(rep(1,N/2), rep(0,N/2))
      u1 <- rep(-0.5, p)
      u2 <- rep( 0.5, p)
      Sigma1 <- matrix(rep(0.5, p^2), nrow = p, ncol = p) + diag(rep(0.5,p))
      Sigma2 <- matrix(rep(0.5, p^2), nrow = p, ncol = p) + diag(rep(0.5,p))
      for(i in 1:p){
        for(j in 1:p){
          Sigma2[i,j] <- 0.5^abs(i-j)
        }
      }
      Xp <- mvrnorm(N/2, u1, Sigma1)
      Xn <- mvrnorm(N/2, u2, Sigma2)
      
      X <- rbind(Xp,Xn)
    }
    X <- cbind(1,X)
  }
  
  Data <- cbind(as.vector(Y), as.matrix(X))
  return(Data)
}