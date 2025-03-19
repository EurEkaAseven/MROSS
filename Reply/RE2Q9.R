rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Reply1")
#################################################
library(MASS)
library(mvtnorm)
library(matrixStats)
library(doParallel) 
library(foreach)

source("DataGenerator.R")
source("Solvers.R")
source("PA.R")
source("getMLE.R")
source("getIBOSS.R")

packs <- c("mvtnorm","MASS","matrixStats","kerndwd")
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
                     method = compare_method[k], wmodel = wmodel, delta = delta,
                     epsilon = C, lambda = lambda, validata = validata)$theta
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
Cs <- 3*log(10)

# Working model
wmodel <- "logistic"

Case <- 1
{# Generating model
  if (Case == 1){
    gmodel <- "logistic"
    option <- "i"
    delta <- 1
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

if (wmodel == "dwd"){
  Cs <-Cs-1 
  compare_method <- c("Unif", "Lopt")
  if(gmodel != "Bdwd"){
    bigDATA <- GenerateData(N = 10*N, p = p, theta_true, option = option, model = gmodel)
    theta_true <- as.vector(getDWD(bigDATA[,-1], bigDATA[,1],
                                   lambda = lambda, validata = validata)$par)}
}else if(gmodel == "X_given_Y"){
  bigDATA <- GenerateData(N = 10*N, p = p, theta_true, option = option, model = gmodel)
  theta_true <- as.vector(getMLE(bigDATA[,-1], bigDATA[,1])$par)}

nm <- length(compare_method)
PA_true <- PreAccuracy(testdata, theta_true)

difn_result <- list()
for (j in 1:length(ns)){
  n <- ns[j]
  system.time({
    cl<- makeCluster(25)  
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
coefEMSE <- EMSE
coefsdMSE <- EMSE

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
    coefMSE <- colSums(bias[seq((k-1)*(p+1)+2,k*(p+1),1),])
    EMSE[j,k] <- mean(MSE, na.rm = T)
    coefEMSE[j,k] <- mean(coefMSE, na.rm = T)
    cat(sum(is.na(MSE)))
    cat(" ")
    sdMSE[j,k] <- sd(MSE, na.rm = T)
    coefsdMSE[j,k] <- sd(coefMSE, na.rm = T)
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

if(wmodel == "logistic" & gmodel != "logistic"){rm("bigDATA")}
save(list = ls(),
     file=paste0("Case",as.character(Case),"_wmodel=",as.character(wmodel),
                 "_gmodel=",as.character(gmodel),
                 "_option=",as.character(option),
                 "_p=", p, "_n0=", n0, "optimal=L_nrep=",nrep,'.RData'))


# plot figure 3 and 4
library(ggplot2)
library(ggsci)
library(cowplot)
# decide which C you want to plot
show.index <- 5
compare_method <- c("Unif", "OSMAC", "IBOSS", "MSCLE")

# intercept
plotMSE <- data.frame(n = rep(ns, nm+1),
                      method = rep(c(compare_method,"MROSS"), each = length(ns)),
                      MSE = log(as.vector(as.matrix(EMSE-coefEMSE)[,c(1:nm,show.index)])),
                      title = paste0("Case ",as.character(Case)))
fig_MSE_intercept <-  
  ggplot(data = plotMSE, aes(x = n, y = MSE, group = method, 
                             color = method, shape = method))+
  geom_point(size=6)+
  geom_line(position = position_dodge(0),cex=1.0)+
  labs(y = "log(MSE)", x = 'r') +
  theme_bw() + 
  theme(
    legend.position = "none",
    #    legend.title = element_text(size = 14),
        legend.title = element_blank(),
       legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) + 
  scale_shape_manual(values = c(8,15,16,17,18)) +
  scale_color_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"))
#  scale_shape_manual(values = c(15,19)) +
#  scale_color_manual(values = c("#ED0000FF","#69C8ECFF"))
dev.new(width=20, height=10, unit="in")
fig_MSE_intercept

# slope
plotMSE <- data.frame(n = rep(ns, nm+1),
                      method = rep(c(compare_method,"MROSS"), each = length(ns)),
                      MSE = log(as.vector(as.matrix(coefEMSE)[,c(1:nm,show.index)])),
                      title = paste0("Case ",as.character(Case)))
fig_MSE_slope <-  
  ggplot(data = plotMSE, aes(x = n, y = MSE, group = method, 
                             color = method, shape = method))+
  geom_point(size=6)+
  geom_line(position = position_dodge(0),cex=1.0)+
  labs(y = "log(MSE)", x = 'r') +
  theme_bw() + 
  theme(
    #legend.position = c(0.85,0.31),
        legend.title = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
  #  legend.background = element_rect(colour ="black"),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)) + 
  scale_shape_manual(values = c(8,15,16,17,18)) +
  scale_color_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"))
#  scale_shape_manual(values = c(15,19)) +
#  scale_color_manual(values = c("#ED0000FF","#69C8ECFF"))
dev.new(width=20, height=10, unit="in")
fig_MSE_slope

re2_fig1 <- plot_grid(fig_MSE_intercept, NULL, fig_MSE_slope, 
          nrow = 1, align = "h", axis = 1, rel_widths = c(1,0.1,1.2), 
          labels=c("A","","B"), label_size = 24)

ggsave(paste0("Intercpt_Case",as.character(Case),"_",as.character(wmodel),"_",
              p,'.png'), 
       re2_fig1, height = 5, width = 14, dpi = 200)