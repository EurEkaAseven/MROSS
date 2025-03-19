rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(matrixStats)
library(doParallel) 
library(foreach)
library(kerndwd)
source("basic/Solvers.R")
source("basic/PA.R")

load("realdata/SUSY_pre.Rdata")
#load("realdata/covtype_pre.Rdata")
#load("realdata/PM25_pre.Rdata")
testdata <- DATA
theta_true <- theta_dwd
packs <- c("mvtnorm","MASS","matrixStats","kerndwd")
simudata_rep <- function(i){
  set.seed(30231102+i)
  ## pilot sampling
  NDATA <- nrow(DATA)
  
  index0 <- sample(NDATA, n0, replace = FALSE)
  # step 1 data
  data0 <- DATA[index0,]
  # step 2 data
  Data <- DATA[-index0,]
  
  all_beta <- vector()
  for (k in 1:nm) {
    set.seed(20231112+i)
    {
        beta <- solver(Data, data0, n, 
                       method = compare_method[k], wmodel = wmodel, 
                       epsilon = C, lambda = lambda, validata = validata)$theta
    }
    if(is.na(beta)){beta <- rep(NA,p+1)}
    all_beta <- c(all_beta,beta)
  }
  for (k in 1:length(Cs)) {
    set.seed(20231112+i)
    beta <- solver(Data, data0, n, 
                   method = "Improve_L", wmodel = wmodel, 
                   Proj_alpha = 1,epsilon = Cs[k], 
                   lambda = lambda, validata = validata)$theta
    if(is.na(beta)){beta <- rep(NA,p+1)}
    all_beta <- c(all_beta,beta)
  }
  
  return(all_beta)
}

nrep <-500                      # Replication size           
n0 <- 1000                      # Step 1 sampling size
ns <- n0*seq(2,5,1)             # Step 2 sampling size
# Working model
wmodel <- "dwd"

  Cs <-seq(1, 4, 1)*log(10)-1
  compare_method <- c("Unif", "Lopt")
  PA_true <- PreAccuracy(testdata, theta_true)
lambda <- 0
validata <- NULL
nm <- length(compare_method)

difn_result <- list()
for (j in 1:length(ns)){
  n <- ns[j]
  system.time({
    cl<- makeCluster(detectCores())  
    registerDoParallel(cl) 
    rep_result <- foreach(
      i=1:nrep,         
      .combine="cbind", 
      #.errorhandling = "pass",
      .packages=packs
    ) %dopar% {simudata_rep(i)}
    stopCluster(cl)
  })
  difn_result[[j]] <- rep_result
}

ntheta <- nm + length(Cs)
EMSE <- matrix(NA,nrow = length(ns),ncol = ntheta)
EMSE <- as.data.frame(EMSE)
rownames(EMSE) <- ns
colnames(EMSE)[1:nm] <- compare_method
colnames(EMSE)[(nm+1):(nm+length(Cs))] <- round(Cs,2)
sdMSE <- EMSE
PA <- EMSE
EL <- EMSE
FS <- EMSE
sdPA <- EMSE
sdEL <- EMSE
sdFS <- EMSE
for (j in 1:length(ns)) {
  bias <- (difn_result[[j]]-theta_true)^2
  for (k in 1:ntheta) {
    MSE <- colSums(bias[seq((k-1)*(p+1)+1,k*(p+1),1),])
    EMSE[j,k] <- mean(MSE, na.rm = T)
    cat(sum(is.na(MSE)))
    cat(" ")
    sdMSE[j,k] <- sd(MSE, na.rm = T)
  }
}

for (j in 1:length(ns)) {
  alltheta <- difn_result[[j]]
  for (k in 1:ntheta) {
    allPA <- colMeans((testdata[,-1]%*%alltheta[seq((k-1)*(p+1)+1,k*(p+1),1),])*(testdata[,1]-0.5) > 0 )
    PA[j,k] <- mean(allPA, na.rm = T)
    sdPA[j,k] <- sd(allPA, na.rm = T)
  }
}
EMSE
PA
PA_true

varBETA <- EMSE
for (j in 1:length(ns)) {
  bias <- (difn_result[[j]]-rowMeans(difn_result[[j]]))^2
  for (k in 1:ntheta) {
    varbeta <- colSums(bias[seq((k-1)*(p+1)+1,k*(p+1),1),])
    varBETA[j,k] <- mean(varbeta, na.rm = TRUE)
  }
}
varBETA


rm("DATA")
rm("testdata")
save(list = ls(),
     file=paste0("SUSY_wmodel=",as.character(wmodel),
                 "_p=", p, "_n0=", n0, "optimal=L",'.RData'))


# plot figures
library(ggplot2)
library(ggsci)
#library(cowplot)
# decide which C you want to plot
show.index <- 5
plotMSE <- data.frame(n = rep(ns, nm+1),
                      method = rep(c(compare_method,"Improve"), each = length(ns)),
                      MSE = log(as.vector(as.matrix(EMSE)[,c(1:nm,show.index)])))
fig_MSE <-  
  ggplot(data = plotMSE, aes(x = n, y = MSE, group = method, 
                             color = method, shape = method))+
  geom_point(size=6)+
  geom_line(position = position_dodge(0),cex=1.0)+
  labs(y = "log(MSE)", x = 'r') +
  theme_bw() + 
  theme(legend.position = "none",
        aspect.ratio = 279/650,
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24)) +
  scale_shape_manual(values = c(15,17,18)) +
  scale_color_manual(values = c("#ED0000FF","#0099B4FF","#925E9FFF"))
dev.new(width=20, height=10, unit="in")
fig_MSE
ggsave(paste0("SUSY_",as.character(wmodel),"_",
              p,'.pdf'), 
       fig_MSE, height = 5, width = 10, dpi = 200)



maxPA <- max(max(PA[,c(1,2,show.index)]),PA_full+0.00005)
minPA <- min(PA[,c(1,2,show.index)])
plotPA <- data.frame(n = rep(ns, nm+1),
                     method = rep(c(compare_method,"Improve"), each = length(ns)),
                     PA = (as.vector(as.matrix(PA)[,c(1:nm,show.index)])))
fig_PA <-  
  ggplot(data = plotPA, aes(x = n, y = PA, group = method, 
                            color = method, shape = method))+
  geom_point(size=6)+
  geom_line(position = position_dodge(0),cex=1.0)+
  labs(y = "Prediction Accuracy", x = 'r') +
  theme_bw() + 
  theme(legend.position = "none",
        aspect.ratio = 279/650,
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24)) +
  scale_shape_manual(values = c(15,17,18)) +
  scale_color_manual(values = c("#ED0000FF","#0099B4FF","#925E9FFF"))+
  geom_abline(slope = 0, intercept = PA_full,linetype = 2,cex = 2) + 
  scale_y_continuous(limits = c(minPA,maxPA),
                     breaks = seq(round(minPA,3)-0.0005,round(maxPA,3),0.0005))

dev.new(width=20, height=10, unit="in")
fig_PA

ggsave(paste0("SUSY_PA_",as.character(wmodel),"_",
              p,'.pdf'), 
       fig_PA, height = 5, width = 10, dpi = 200)