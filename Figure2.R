rm(list = ls())
gc()
######## set your working directory here ########
setwd("D:/BIT/MASTER/code/Codes_Reproducibility")
#################################################library(mvtnorm)
library(MASS)
library(ggplot2)
library(OptimalDesign)
library(ggsci)
library(directlabels)
library(latex2exp)
source("basic/DataGenerator.R")
source("basic/Solvers.R")
source("basic/getMLE.R")

# basic settings
N <- 1e5                       # Full data size
p <- 2                         # Dimension of X
n <- 500                       # sampling size
n0 <- 500                       # pilot sampling size

wmodel <- "logistic"
gmodel <- "logistic"
option <- "vii"
radius_unif <- 6
plot_range <- 6

# gennerating full sample
set.seed(20240328)
theta_true <- c(0,rep(0.5,p))      # True theta of logistic generating model
data0 <- GenerateData(N = n0, p = p, theta_true, 
                      option = option, model = gmodel)
theta_pilot <- as.vector(getMLE(data0[,-1], data0[,1])$par)
DATA <- GenerateData(N = N, p = p, theta_true, 
                     option = option, model = gmodel)

full_point <- data.frame(X1 = DATA[,3], X2 = DATA[,4], 
                         Y = as.factor(2*DATA[,1]-1))

# pilot data
dataX <- DATA[,-1]
dataY <- DATA[, 1]
dataN <- nrow(DATA)
C <- 2.5
#pilot_theta <- as.vector(getMLE(data0[,-1], data0[,1])$par)
sample_result <- solver(DATA, data0, n, method = "Improve_A", wmodel = "logistic",
                        epsilon = C)

centroid <- sample_result$xbar
centroid <- data.frame(cY = as.factor(c(1,0)), cX1 = centroid[2,], cX2 = centroid[3,])
sample_point <- sample_result$sdata
weight_OSMAC <- (sample_result$weight)/(sample_result$weight_p)
weight <-  sample_result$weight
weight_p <- sample_result$weight_p
sample_point <- data.frame(sY = sample_point[,1], 
                           sX1 = sample_point[,3],
                           sX2 = sample_point[,4])


# confidence ellipse
# mu is mean, s is covariance matrix, c is a given number associated with confidence level
ellipse.simple <- function(s1, s2, c){
  a <- s1*c
  b <- s2*c 
  x <- seq(from = -a, to = a, length.out = 400)
  points <- data.frame(
    x1 = c(-x, x),
    x2 = NA
  )
  points$x2[1:400] <- sqrt(((a*b)^2-(b*x)^2)/a^2)
  points$x2[401:800] <- -sqrt(((a*b)^2-(b*x)^2)/a^2)
  return(points)
}
ellipse.general <- function(mu, s, c){
  lambda <- diag(eigen(s)$values)  # eigen value
  P <- eigen(s)$vectors            # eigen vector
  Y <- ellipse.simple(s1 = sqrt(lambda[1,1]), s2 = sqrt(lambda[2,2]), c = c)  # 中心在原点，没有倾斜的椭圆的坐标
  X <- t(P%*%t(Y) + mu)
  X <- as.data.frame(X)
  colnames(X) <- c('x1', 'x2')
  return(X)
}

mu <- 0  
s <- matrix(rep(0.5, p^2), nrow = p, ncol = p) + diag(rep(0.5,p))
s <- sqrt(2)*s
c_alpha <- 0.001
c <- sqrt(2*log(1/(c_alpha)))

X_ellipse <- ellipse.general(mu = mu, s = s, c = c)
X_ellipse <- data.frame(X_e1 = X_ellipse[,1], X_e2 = X_ellipse[,2])


form.quad <- ~ x1 + x2 
Fx <- Fx_glm(form.quad, theta_true,
             glm.model="bin-logit", lower = c(-radius_unif,-radius_unif), 
             upper = c(radius_unif,radius_unif), n.levels=c(100,100))

Fx.lin <- Fx_cube(form.quad , lower = c(-radius_unif,-radius_unif), 
                  upper = c(radius_unif,radius_unif), n.levels=c(100,100)) # Just for the plot

# design points
lambda <- eigen(s)$values
a <- lambda[1]*c^2
b <- lambda[2]*c^2
Trans <- eigen(s)$vectors
tFx.lin <- t(t(Trans)%*%t(Fx.lin[,2:3])) + mu
design_index <- which(tFx.lin[, 1]^2/a + tFx.lin[, 2]^2/b - 1 < 0)

design_point <- Fx.lin[design_index, 2:3]
Fx <- Fx[design_index,]
w <- od_REX(Fx, crit = "A")$w.best

#od_plot(Fx, w, design_point, dd.size=2)
suppot_point <- design_point[w>0.1,]
suppot_point <- data.frame(supX1 = suppot_point[,1],
                           supX2 = suppot_point[,2])
dd <- dirder(Fx, w)
design_point <- data.frame(dX1 = design_point[,1],
                           dX2 = design_point[,2],
                           dd = dirder(Fx, w))

#breaks_lines <- round(quantile(dd, seq(0.1, 0.9, 0.2)), 1)
breaks_lines <- seq(-0.6,-0.2, 0.4)

# layer of suppot point of optimal design
layer_suppot <- geom_point(data = suppot_point,aes(x = supX1, y = supX2),
                           color = "red",shape = 4, cex = 5,stroke = 2)

# layer of confidence ellipse
layer_ellipse <- geom_path(data = X_ellipse, 
                           mapping = aes(x = X_e1, y = X_e2), 
                           color = 'black', linetype = 6)

fig_design <- ggplot(data = design_point, 
                     aes(x = dX1, y = dX2)) +
  xlim(-plot_range,plot_range) + ylim(-plot_range,plot_range) +
  labs(x = expression(X[1]), y = expression(X[2])) +
  theme_classic() + 
  coord_fixed()


# Fig2(a)
fig1 <- ggplot() +
  geom_point(data = full_point, aes(x = X1, y = X2), shape = 1,
             color = "grey80", size = 1) +
  xlim(-plot_range,plot_range) + ylim(-plot_range,plot_range) +
  labs(x = expression(X[1]), y = expression(X[2]), 
       shape = "Full sample points", color = "dd") +
  theme_classic() + 
  coord_fixed() + 
  stat_contour(data = design_point, aes(x = dX1, y = dX2, z = dd, 
                                        colour = ..level..), 
               breaks=breaks_lines, cex = 0.8) +
  scale_color_gradient(low = "grey20", high = "grey20") +
  layer_ellipse +
  layer_suppot +
  geom_point(data = centroid[1,], aes(x = cX1, y = cX2),
             size = 5, color = "red", shape = 17) +
  geom_point(data = centroid[2,], aes(x = cX1, y = cX2),
             size = 5, color = "red", shape = 16) + 
  annotate("text", x= 3 , y= 5.5 , label= expression(chi["+"]),
           size = 14, family = "serif") +
  annotate("text", x= -5 , y= -4 , label= expression(chi["-"]),
           size = 14, family = "serif") +
  annotate("text", x= -3 , y= 5 , label= expression(chi["s"]),
           size = 14, family = "serif") +
  geom_abline(slope = -theta_pilot[3]/theta_pilot[2],
              intercept = -(theta_pilot[1]+C)/theta_pilot[2],
              linetype = 5, cex = 1) +
  geom_abline(slope = -theta_pilot[3]/theta_pilot[2],
              intercept = -(theta_pilot[1]-C)/theta_pilot[2],
              linetype = 5, cex = 1) + 
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -theta_true[1]/theta_true[2],
              linetype = 1, cex = 1.2) +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) 


dev.new(width=20, height=20, unit="in")
direct.label(fig1, list("bottom.pieces", colour='red', 
                        cex = 2,
                        fontface="plain", fontfamily="serif"))




plot_DATA <-  sample_result$tdata
expcomp <- exp(plot_DATA[,-1]%*%theta_pilot)
P <- 1-1/(1+expcomp)
H <- t(plot_DATA[,-1])%*%(as.vector(P * (1-P)) * plot_DATA[,-1])
Pi_tmp <- abs(plot_DATA[,1]-P) * sqrt(colSums(solve(H, t(plot_DATA[,-1]))^2))
#Pi_tmp <- abs(plot_DATA[,1]-P) * sqrt(rowSums(plot_DATA[,-1]^2))
Pi <- log(pmin(1,n*Pi_tmp/sum(Pi_tmp)))
plot_DATA <- data.frame(pX1 = plot_DATA[,3],
                        pX2 = plot_DATA[,4],
                        Pi = Pi)
# bining
radius_unif <- 6
delta <- 0.3
design_point <- expand.grid(X = seq(-radius_unif ,radius_unif ,delta),
                            Y = seq(-radius_unif ,radius_unif ,delta))
dim_lattice <- 2*radius_unif/delta
Pi_sum <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
Pi_num <- matrix(0, nrow = dim_lattice+1, ncol = dim_lattice+1)
grid1 <- floor(plot_DATA$pX1/delta) + radius_unif/delta +1
grid2 <- floor(plot_DATA$pX2/delta) + radius_unif/delta +1
for (i in 1:nrow(plot_DATA)) {
  Pi_sum[grid1[i], grid2[i]] <-  Pi_sum[grid1[i], grid2[i]] + Pi[i]
  Pi_num[grid1[i], grid2[i]] <- Pi_num[grid1[i], grid2[i]] + 1
}

Pi_mean <- Pi_sum/Pi_num

Pi_mean <- as.vector(Pi_mean)
low_index <- which(as.vector(Pi_num)<4)

Pi_mean <- Pi_mean[-low_index]
design_point <- design_point[-low_index,]

design_point$Z <- Pi_mean
design_point$X <- design_point$X + delta/2
design_point$Y <- design_point$Y + delta/2
fig_heatmap <- ggplot(design_point,aes(x=X,y=Y))+
  geom_tile(aes(fill = Z), colour = "white") +
  scale_fill_material("grey") +
  xlim(-plot_range,plot_range) + ylim(-plot_range,plot_range) +
  theme_classic() +
  coord_fixed()

# Fig2(b)
fig2 <- fig_heatmap +
  layer_suppot + 
  layer_ellipse + 
  geom_point(data = centroid[1,], aes(x = cX1, y = cX2),
             size = 5, color = "red", shape = 17) +
  geom_point(data = centroid[2,], aes(x = cX1, y = cX2),
             size = 5, color = "red", shape = 16) +
  annotate("text", x= 3 , y= 5.5 , label= expression(chi["+"]),
           size = 14, family = "serif") +
  annotate("text", x= -5 , y= -4 , label= expression(chi["-"]),
           size = 14, family = "serif") +
  annotate("text", x= -3 , y= 5 , label= expression(chi["s"]),
           size = 14, family = "serif") +
  geom_abline(slope = -theta_pilot[3]/theta_pilot[2],
              intercept = -(theta_pilot[1]+C)/theta_pilot[2],
              linetype = 5, cex = 1) +
  geom_abline(slope = -theta_pilot[3]/theta_pilot[2],
              intercept = -(theta_pilot[1]-C)/theta_pilot[2],
              linetype = 5, cex = 1) + 
  geom_abline(slope = -theta_true[3]/theta_true[2],
              intercept = -theta_true[1]/theta_true[2],
              linetype = 1, cex = 1.2) +
  theme(
    legend.position = c(0.95,0.85),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  labs(x = expression(X[1]), y = expression(X[2]),
       fill = TeX(r'($\log(\tilde{\pi}^{A})$)') )

dev.new(width=30, height=20, unit="in")
fig2


# save the figure
ggsave("normal_design.png", 
       direct.label(fig1, list("bottom.pieces", colour='red', cex=2, 
                               fontface="plain", fontfamily="serif")) ,
       height = 7, width = 8, dpi = 200)
ggsave("normal_design_prob.png", fig2, height = 7, width = 8, dpi = 200)