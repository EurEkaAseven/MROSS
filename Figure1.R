rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(MASS)
library(mvtnorm)
library(ggplot2)
library(kerndwd)
library(ggsci)
library(latex2exp)
source("basic/DataGenerator.R")
source("basic/getMLE.R")
source("basic/getsSVM.R")

set.seed(20240328)      
N <- 1e5                      # Full data size
p <- 2                        # Dimension of X
n0 <- 200                     # Step 1 sampling size
ns <- 500                     # Step 2 sampling size

# Generating model 
gmodel <- "logistic" 
option <- "vii"               # Distribution setting of X
theta_true <- c(0,rep(0.5,p)) # True theta of logistic generating model

# tuning parameters for DWD and sSVM
validata <- NULL
lambda <- 0
Cost <- 1e-4

Data <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
X <- Data[,-1]
Y <- Data[, 1]
Xnorm <- sqrt(rowSums(X^2))
theta_logistic <- as.vector(getMLE(X, Y)$par)
theta_dwd <- as.vector(kerndwd(X[,-1], Y, lambda = lambda, kern = vanilladot())$alpha)
theta_ssvm <- as.vector(getsSVM(X, Y, C=Cost)$par)

############## Fig1(a) ################
C <- 2.7
u <- X%*%theta_logistic
P <- 1-1/(1+exp(u))
g <- abs(Y-P)
Prob <- g*Xnorm
Prob <- log(pmin(1,ns*Prob/sum(Prob)))
theta_true <- theta_logistic
# binning
radius_unif <- 6
delta <- 0.3
design_point <- expand.grid(X = seq(-radius_unif ,radius_unif ,delta),
                            Y = seq(-radius_unif ,radius_unif ,delta))
dim_lattice <- 2*radius_unif/delta
Pi_sum <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
Pi_num <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
grid1 <- floor(X[,2]/delta) + radius_unif/delta +1
grid2 <- floor(X[,3]/delta) + radius_unif/delta +1
for (i in 1:nrow(X)) {
  Pi_sum[grid1[i], grid2[i]] <-  Pi_sum[grid1[i], grid2[i]] + Prob[i]
  Pi_num[grid1[i], grid2[i]] <- Pi_num[grid1[i], grid2[i]] + 1
}

Pi_mean <- Pi_sum/Pi_num
Pi_mean <- as.vector(Pi_mean)
low_index <- which(as.vector(Pi_num)<3)

Pi_mean <- Pi_mean[-low_index]
design_point <- design_point[-low_index,]

design_point$Z <- Pi_mean
design_point$X <- design_point$X + delta/2
design_point$Y <- design_point$Y + delta/2

logistic_bin_prob <- ggplot(design_point,aes(x=X,y=Y))+
  geom_tile(aes(fill = Z), colour = "white") +
  scale_fill_material("grey") +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -theta_true[1]/theta_true[2],
              linetype = 1, cex = 1.2)+
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]+C)/theta_true[2],
              linetype = 5, cex = 1) +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]-C)/theta_true[2],
              linetype = 5, cex = 1) +
  annotate("text", x= 5 , y= 5 , label= expression(chi["+"]),
           size = 14, family = "serif") +
  annotate("text", x= -5 , y= -5 , label= expression(chi["-"]),
           size = 14, family = "serif") +
  annotate("text", x= -4 , y= 5 , label= expression(chi["s"]),
           size = 14, family = "serif") +
  xlim(-6,6) + ylim(-6,6) +
  theme_classic() +
  coord_fixed() + 
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(x = expression(X[1]), y = expression(X[2]),
       fill = TeX(r'($\log(\pi^{L})$)'))
dev.new(width=20, height=20, unit="in")
logistic_bin_prob
#######################################

############## Fig1(b) ################
C <- 2.3
u <- (X%*%theta_dwd)*(2*Y-1)
g <- ifelse(u<=0.5, 1, 1/u^2/4)
Prob <- g*Xnorm
Prob <- log(pmin(1,ns*Prob/sum(Prob)))
theta_true <- theta_dwd

# binning
radius_unif <- 6
design_point <- expand.grid(X = seq(-radius_unif ,radius_unif ,delta),
                            Y = seq(-radius_unif ,radius_unif ,delta))
dim_lattice <- 2*radius_unif/delta
Pi_sum <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
Pi_num <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
grid1 <- floor(X[,2]/delta) + radius_unif/delta +1
grid2 <- floor(X[,3]/delta) + radius_unif/delta +1
for (i in 1:nrow(X)) {
  Pi_sum[grid1[i], grid2[i]] <-  Pi_sum[grid1[i], grid2[i]] + Prob[i]
  Pi_num[grid1[i], grid2[i]] <- Pi_num[grid1[i], grid2[i]] + 1
}

Pi_mean <- Pi_sum/Pi_num

Pi_mean <- as.vector(Pi_mean)
low_index <- which(as.vector(Pi_num)<3)

Pi_mean <- Pi_mean[-low_index]
design_point <- design_point[-low_index,]

design_point$Z <- Pi_mean
design_point$X <- design_point$X + delta/2
design_point$Y <- design_point$Y + delta/2

dwd_bin_prob <- ggplot(design_point,aes(x=X,y=Y))+
  geom_tile(aes(fill = Z), colour = "white") +
  scale_fill_material("grey") +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -theta_true[1]/theta_true[2],
              linetype = 1, cex = 1.2)+
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]+C)/theta_true[2],
              linetype = 5, cex = 1) +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]-C)/theta_true[2],
              linetype = 5, cex = 1) +
  annotate("text", x= 5 , y= 5 , label= expression(chi["+"]),
           size = 14, family = "serif") +
  annotate("text", x= -5 , y= -5 , label= expression(chi["-"]),
           size = 14, family = "serif") +
  annotate("text", x= -4 , y= 5 , label= expression(chi["s"]),
           size = 14, family = "serif") +
  xlim(-6,6) + ylim(-6,6) +
  theme_classic() +
  coord_fixed() + 
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(x = expression(X[1]), y = expression(X[2]),
       fill = TeX(r'($\log(\pi^{L})$)'))
dev.new(width=20, height=20, unit="in")
dwd_bin_prob
#######################################

############### Fig1(c) ###############
C <- 1.2
u <- (X%*%theta_ssvm)*(2*Y-1)
g <- ifelse(1-u > 0, 2*Cost*(1-u), 0)
Prob <- g*Xnorm + 0.0001
Prob <- log(pmin(1,ns*Prob/sum(Prob)))
theta_true <- theta_ssvm
# binning
radius_unif <- 6
design_point <- expand.grid(X = seq(-radius_unif ,radius_unif ,delta),
                            Y = seq(-radius_unif ,radius_unif ,delta))
dim_lattice <- 2*radius_unif/delta
Pi_sum <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
Pi_num <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
grid1 <- floor(X[,2]/delta) + radius_unif/delta +1
grid2 <- floor(X[,3]/delta) + radius_unif/delta +1
for (i in 1:nrow(X)) {
  Pi_sum[grid1[i], grid2[i]] <-  Pi_sum[grid1[i], grid2[i]] + Prob[i]
  Pi_num[grid1[i], grid2[i]] <- Pi_num[grid1[i], grid2[i]] + 1
}

Pi_mean <- Pi_sum/Pi_num

Pi_mean <- as.vector(Pi_mean)
low_index <- which(as.vector(Pi_num)<3)

Pi_mean <- Pi_mean[-low_index]
design_point <- design_point[-low_index,]

design_point$Z <- Pi_mean
design_point$X <- design_point$X + delta/2
design_point$Y <- design_point$Y + delta/2
ssvm_bin_prob <- ggplot(design_point,aes(x=X,y=Y))+
  geom_tile(aes(fill = Z), colour = "white") +
  scale_fill_material("grey") +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -theta_true[1]/theta_true[2],
              linetype = 1, cex = 1.2)+
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]+C)/theta_true[2],
              linetype = 5, cex = 1) +
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -(theta_true[1]-C)/theta_true[2],
              linetype = 5, cex = 1) +
  annotate("text", x= 5 , y= 5 , label= expression(chi["+"]),
           size = 14, family = "serif") +
  annotate("text", x= -5 , y= -5 , label= expression(chi["-"]),
           size = 14, family = "serif") +
  annotate("text", x= -4 , y= 5 , label= expression(chi["s"]),
           size = 14, family = "serif") +
  xlim(-6,6) + ylim(-6,6) +
  theme_classic() +
  coord_fixed() + 
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))  +
  labs(x = expression(X[1]), y = expression(X[2]),
       fill = TeX(r'($\log(\pi^{L})$)'))
dev.new(width=20, height=20, unit="in")
ssvm_bin_prob
######################################

# save the figures
ggsave("logistic.png", logistic_bin_prob , width = 7, height = 7, dpi = 200)
ggsave("dwd.png",dwd_bin_prob , width = 7, height = 7, dpi = 200)
ggsave("ssvm.png", ssvm_bin_prob , width = 7, height = 7, dpi = 200)