#2019/11/02

#comparison of predictive performances by using simulation data

rm(list=ls(all = T));
gc();gc();           

##################################################################################################
###     package
##################################################################################################

#=================================================================================================
###   for Monte Carlo simulation
#=================================================================================================

library(magrittr)  #pipe
library(tcltk)     #ProgressBar

#=================================================================================================
###   for plot
#=================================================================================================

library(ggplot2)
library(reshape2)

##################################################################################################
###     setting
##################################################################################################

n <- 50                  #sample size
ks <- c((5:9)*n/10, n-4)  #the numbers of explanatory variables
rho <- 0.99               #parameter of correlation between explanatory variables
ms.n <- 100             #iteration of Monte Carlo simulation
case <- 1                 #change definition of beta

##################################################################################################
###     path
##################################################################################################

# setwd()

##################################################################################################
###     function
##################################################################################################

#=================================================================================================
###   GCp minimization method
#=================================================================================================

source("f_GRRGCp.R")

#=================================================================================================
###   other
#=================================================================================================

# calculating square root of symmetric matrix

#argument
#  M: symmetrix matrix

mat.root <- function(M)
{
  k <- nrow(M)
  M.eig <- eigen(M)
  Q <- M.eig$vectors
  D <- M.eig$values %>% sqrt %>% diag(., k, k)
  Out <- Q %*% D %*% t(Q)
  return(Out)
}

##################################################################################################
###     main
##################################################################################################

RMSE <- TIME <- matrix(0, length(ks), 5)
ALPHA <- matrix(0, ms.n, length(ks))

pb <- txtProgressBar(min = 1, max = length(ks)*ms.n, style = 3)

for(k.i in 1:length(ks))
{
  k <- ks[k.i]
  
  ##################################################################################################
  ###     create data
  ##################################################################################################
  
  #=================================================================================================
  ###   X
  #=================================================================================================
  
  set.seed(123)  #fix random variables
  X0 <- runif(n*k, -1, 1) %>% matrix(., n, k)
  
  In <- diag(n)
  Jn <- matrix(1/n, n, n)
  
  PHI <- matrix(0, k, k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      PHI[i,j] <- PHI[j,i] <- rho^abs(i-j)
    }
  }
  PHI2 <- mat.root(PHI)
  
  X <- (In - Jn) %*% X0 %*% PHI2
  
  #=================================================================================================
  ###   beta
  #=================================================================================================
  
  if(case == 1)
  {
    Beta <- rep(1, k)
  } else if(case == 2)
  {
    Beta <- (1:k)/k
  } else
  {
    stop("this case is not defined")
  }
  
  
  ##################################################################################################
  ###     other preparation
  ##################################################################################################
  
  XBeta <- X %*% Beta
  
  t.s0 <- proc.time()[3]
  X.svd <- svd(X)
  t.f0 <- proc.time()[3]
  
  SE <- matrix(0, ms.n, 5)
  Time <- numeric(5)
  
  ##################################################################################################
  ###     Monte Carlo Simulation
  ##################################################################################################
  
  for(ms.i in 1:ms.n)
  {
    setTxtProgressBar(pb, (k.i-1)*ms.n + ms.i)
    
    t.s <- t.f <- numeric(5)
    
    set.seed(ms.i)
    y <- XBeta + rnorm(n)
    
    t.s00 <- proc.time()[3]
    res <- GRRGCp(y, X, centralize=T, X.svd=X.svd)
    t.f00 <- proc.time()[3]
    
    ##################################################################################################
    ###     Optimal Cp
    ##################################################################################################
    
    t.s[1] <- proc.time()[3]
    est1 <- get.est(res)
    t.f[1] <- proc.time()[3]
    
    ALPHA[ms.i, k.i] <- est1$alpha
    yh1 <- est1$y.hat
    SE[ms.i, 1] <- sum((yh1 - XBeta)^2)
    
    ##################################################################################################
    ###     alpha = 1 (Cp)
    ##################################################################################################

    t.s[2] <- proc.time()[3]
    est2 <- get.est(res, alpha=1)
    t.f[2] <- proc.time()[3]
    
    yh2 <- est2$y.hat
    SE[ms.i, 2] <- sum((yh2 - XBeta)^2)
    
    ##################################################################################################
    ###     alpha = 1 + (2/(n-k-3)) (MCp)
    ##################################################################################################
    
    t.s[3] <- proc.time()[3]
    est3 <- get.est(res, alpha=1 + (2/(n-k-3)))
    t.f[3] <- proc.time()[3]
    
    yh3 <- est3$y.hat
    SE[ms.i, 3] <- sum((yh3 - XBeta)^2)
    
    ##################################################################################################
    ###     alpha = log log n (HQC-type Cp)
    ##################################################################################################

    t.s[4] <- proc.time()[3]
    est4 <- get.est(res, alpha=log(log(n)))
    t.f[4] <- proc.time()[3]
    
    yh4 <- est4$y.hat
    SE[ms.i, 4] <- sum((yh4 - XBeta)^2)
    
    ##################################################################################################
    ###     alpha = log n / 2 (BIC-type Cp)
    ##################################################################################################

    t.s[5] <- proc.time()[3]
    est5 <- get.est(res, alpha=log(n)/2)
    t.f[5] <- proc.time()[3]
    
    yh5 <- est5$y.hat
    SE[ms.i, 5] <- sum((yh5 - XBeta)^2)
    
    ##################################################################################################
    ###     time
    ##################################################################################################
    
    Time <- Time + (t.f00 - t.s00) + (t.f - t.s)
    
  } #end for ms.i
  
  MSE <- apply(SE, 2, mean)
  RMSE[k.i, ] <- 100 * MSE/(k+1)
  TIME[k.i, ] <- (t.f0 - t.s0) + (Time/ms.n)
} #end for k.i

##################################################################################################
###     summary
##################################################################################################

RMSE <- data.frame(RMSE)
TIME <- data.frame(TIME)
colnames(RMSE) <- colnames(TIME) <- c("Optimal", "Cp", "MCp", "HQC-type", "BIC-type")
rownames(RMSE) <- rownames(TIME) <- ks

colnames(ALPHA) <- ks

ALPHAs <- cbind(apply(ALPHA, 2, mean), 1, 1 + (2/(n-ks-3)), log(log(n)), log(n)/2)
colnames(ALPHAs) <- c("Optimal", "Cp", "MCp", "HQC-type", "BIC-type")

##################################################################################################
###     figure
##################################################################################################

#=================================================================================================
###  MSE
#=================================================================================================

data <- data.matrix(RMSE)
data <- melt(data)
colnames(data) <- c("k", "criterion", "RMSE")

fig1 <- ggplot(data=data, aes(x=k, y=RMSE, shape=criterion)) +
  geom_line(size=1) + 
  geom_point(size=5) + 
  xlab("the number of explanatory variables") + 
  ylab("Relative MSE") +
  theme(
    axis.title.x = element_text(size=20), 
    axis.title.y = element_text(size=20), 
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20),
    legend.title=element_blank(),
    legend.text = element_text(size=20),
    legend.position = "top",
    legend.justification = c(0,0),
    plot.margin= unit(c(0, 1, 0, 0), "lines")
  ) 

#=================================================================================================
###  alpha
#=================================================================================================

data <- melt(ALPHA)
colnames(data) <- c("No", "k", "alpha")
data$k <- data$k %>% as.character

fig2 <- ggplot(data, aes(x = k, y = alpha)) +
  geom_boxplot() +
  xlab("the number of explanatory variables") + 
  ylab("optimal alpha")+
  theme(
    axis.title.x = element_text(size=18), 
    axis.title.y = element_text(size=18), 
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18)
  ) 

#=================================================================================================
###  running time
#=================================================================================================

data <- data.matrix(TIME)
data <- melt(data)
colnames(data) <- c("k", "criterion", "time")

fig3 <- ggplot(data=data, aes(x=k, y=time, shape=criterion)) +
  geom_line(size=1) + 
  geom_point(size=4) + 
  xlab("the number of explanatory variables") + 
  ylab("running time (s)") +
  theme(
    axis.title.x = element_text(size=18), 
    axis.title.y = element_text(size=18), 
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18)
  )  