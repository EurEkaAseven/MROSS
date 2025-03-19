rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################
library(MASS)
library(mvtnorm)
library(matrixStats)
library(doParallel) 
library(foreach)
source("basic/DataGenerator.R")
source("basic/Solvers.R")
source("basic/PA.R")
source("basic/getIBOSS.R")

packs <- c("mvtnorm","MASS","matrixStats")
simudata_rep <- function(i){
  set.seed(30231102+i)
  DATA <- GenerateData(N = N, p = p, theta_true, option = option, model = gmodel)
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
                     epsilon = C, lambda = lambda, validata = validata, 
                     delta = delta)$theta
    }
    if(is.na(beta)){beta <- rep(NA,p+1)}
    all_beta <- c(all_beta,beta)
  }
  for (k in 1:length(Cs)) {
    set.seed(20231112+i)
    beta <- solver(Data, data0, n, 
                   method = "Improve_L", wmodel = wmodel, 
                   epsilon = Cs[k], 
                   lambda = lambda, validata = validata)$theta
    if(is.na(beta)){beta <- rep(NA,p+1)}
    all_beta <- c(all_beta,beta)
  }
  return(all_beta)
}
# basic settings
nrep <-500                    # Replication number
N <- 5e5                      # Full data size
p <- 20                       # Dimension of X (without intercept term)
n0 <- 1000                    # Step 1 sampling size
ns <- n0*seq(2,5,1)           # Step 2 sampling size
Cs <- seq(1, 4, 1)*log(10)    # varying Cr0

# Working model
wmodel <- "logistic"

Case <- 1
{# Generating model
  if (Case == 1){
    gmodel <- "logistic"
    option <- "i"
    delta <- 1    # the tuning parameter of IBOSS
  }else if(Case == 2){
    gmodel <- "logistic"
    option <- "ii"
    delta <- 1.5
  }else if(Case == 3){
    gmodel <- "X_given_Y"
    option <- "vi"
    delta <- 10
  }else if(Case == 4){
    gmodel <- "X_given_Y"
    option <- "v"
    delta <- 10
  }
}

set.seed(20231102)
theta_true <- c(0.5,rep(0.5,p))  # True theta of logistic generating model

# Test set
testdata <- GenerateData(N = 0.1*N, p = p, theta_true, 
                         option = option, model = gmodel)
# parameters for dwd
validata <- NULL
lambda <- 0

compare_method <- c("Unif", "OSMAC_L", "IBOSS", "MoreEff_Wang_L")
nm <- length(compare_method)
PA_true <- PreAccuracy(testdata, theta_true)

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

if(gmodel != "logistic"){
  bigDATA <- GenerateData(N = 10*N, p = p, theta_true, option = option, model = gmodel)
  theta_true <- as.vector(getDWD(bigDATA[,-1], bigDATA[,1],
                                 lambda = lambda, validata = validata)$par)}

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

varBETA <- EMSE
for (j in 1:length(ns)) {
  bias <- (difn_result[[j]]-rowMeans(difn_result[[j]]))^2
  for (k in 1:ntheta) {
    varbeta <- colSums(bias[seq((k-1)*(p+1)+1,k*(p+1),1),])
    varBETA[j,k] <- mean(varbeta, na.rm = TRUE)
  }
}


if(wmodel == "logistic" & gmodel != "logistic"){rm("bigDATA")}
save(list = ls(),
     file=paste0("Case",as,character(Case),"_wmodel=",as.character(wmodel),
                 "_gmodel=",as.character(gmodel),
                 "_option=",as.character(option),
                 "_p=", p, "_n0=", n0, "optimal=L_nrep=",nrep,'.RData'))

# plot figure 3 and 4
library(ggplot2)
library(ggsci)

show.index <- 6  # decide which Cr0 you want to plot
compare_method <- c("AUnif", "BOSMAC_L", "DIBOSS", "CMoreEff")
# Figure 3(a)
plotMSE <- data.frame(n = rep(ns, nm+1),
                      method = rep(c(compare_method,"EImprove"), each = length(ns)),
                      MSE = log(as.vector(as.matrix(EMSE)[,c(1:nm,show.index)])))
fig_MSE <-  
  ggplot(data = plotMSE, aes(x = n, y = MSE, group = method, 
                             color = method, shape = method))+
  geom_point(size=6)+
  geom_line(position = position_dodge(0),cex=1.0)+
  labs(y = "log(MSE)", x = 'r') +
  theme_bw() + 
  theme(legend.position = "none",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) + 
  scale_shape_manual(values = c(18,17,16,8,15)) +
  scale_color_manual(values = c("#925E9FFF","#0099B4FF","#42B540FF","#00468BFF","#ED0000FF"))
dev.new(width=20, height=10, unit="in")
fig_MSE
ggsave(paste0("Case",as.character(Case),"_",as.character(wmodel),"_",
              p,'.pdf'), 
       fig_MSE, height = 5, width = 10, dpi = 200)

# Figure 4(a)
mse <- vector()
for (j in c(1,length(ns))) {
  bias <- (difn_result[[j]]-theta_true)^2
  for (k in 1:ntheta) {
    if( (k <= nm) | k == show.index){
      mse <- c(mse,colSums(bias[seq((k-1)*(p+1)+1,k*(p+1),1),]))
    }
  }
}

boxplot_mse <- data.frame(r = as.factor(rep(c(ns[1],ns[length(ns)]),each = nrep*(nm+1))),
                          method = rep(rep(c(compare_method, "EImprove"),each = nrep),2),
                          MSE = mse)

fig_box <- ggplot(boxplot_mse,aes(x = r,y = MSE, color = method))+
  stat_boxplot(geom="errorbar",position=position_dodge(0.6),width = 0.5)+
  geom_boxplot(position=position_dodge(0.6),
               size=0.2,
               width=0.6,
               outlier.shape = NA,
  )+
  scale_color_manual(values = c("#925E9FFF","#0099B4FF","#42B540FF","#00468BFF","#ED0000FF"))+
ylim(0, 0.25) +
  theme_classic() + 
  theme(
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20))
dev.new(width=20, height=10, unit="in")
fig_box
ggsave(paste0("Case",as.character(Case),"box_",as.character(wmodel),"_",
              p,'.pdf'), 
       fig_box, height = 7, width = 7, dpi = 200)