#2019/11/03

#These functions calculate generalized ridge regression estimator based on GCp criterion minimimzation method.

#####################################################################################################

#This function gives what required for obtaining estimation result

#argumant
#  y: n-dimensional vector of a response variable
#  X: n \times k matrix of an explanatory variable; k < n - 3
#  centralize: whether X is centralized or not
#  X.svd: result of "svd(X)"; default is FALSE

#output
#  what required for obtaining estimation results

GRRGCp <- function(y, X, centralize=F, X.svd=F)
{
  n <- length(y)
  k <- ncol(X)
  q <- 1 - (2/(n-k-1))
  
  if(q <= 0){stop("q is nonpositive")}
  if(n-1 <= k){stop("GCp cannot be defined")}
  
  if(centralize == F)
  {
    In <- diag(n)
    Jn <- matrix(1/n, n, n)
    X <- (In - Jn) %*% X
  }
  
  if(class(X.svd) == "logical"){X.svd <- svd(X)}
  P1 <- X.svd$u
  P1. <- t(P1)
  d2 <- X.svd$d
  # if(identical(all.equal(d2[k], 0), T)){stop("X is rank deficient")}
  if(d2[k] < 1e-7){stop("X is rank deficient")}
  Q <- X.svd$v
  
  z <- as.vector(P1. %*% y)
  z2 <- z^2
  
  s0 <- ( sum(y^2) - (sum(y)^2)/n - sum(z2) ) / (n-k-1)
  
  Q.Beta <- (1/d2)*z
  
  mu <- mean(y)
  s.t <- sum((y - mu)^2)
  
  Out <- list(
    n = n,
    k = k,
    q = q,
    z = z,
    z2 = z2,
    s0 = s0,
    d2 = d2,
    Q = Q,
    Q.Beta = Q.Beta,
    mu = mu,
    X = X,
    y = y,
    s.t = s.t
  )
  
  return(Out)
}

#####################################################################################################

#This function gives estimation results at given alpha via "GRRGCp".

#argument
#  res: output of GRRGCp()
#  alpha: penalty in GCp criterion; numeric or "optimal" (default) ; "optimal" uses alpha optimized by the proposed method

#output
#  y.hat: vector of predictive values
#  beta.hat: estimator of beta
#  R2

get.est <- function(res, alpha="optimal")
{
  k <- res$k
  q <- res$q
  z2 <- res$z2
  s0 <- res$s0
  Q <- res$Q
  Q.Beta <- res$Q.Beta
  mu <- res$mu
  X <- res$X
  y <- res$y
  s.t <- res$s.t
  
  if(alpha == "optimal")
  {
    t0 <- z2/s0
    t <- c(0, sort(t0))
    c1 <- cumsum(t)
    c2 <- rev(cumsum(c(0, rev(1/t[-1]))))
    
    Phi <- q*c2*(t^2) + 2*c2*t - 2*(0:k) + q*c1
    
    alpha <- t[which.min(Phi)]
  }
  
  v <- as.vector(1 - ((alpha*s0)/(z2)))
  v[v < 0] <- 0
  
  Betah <- Q %*% (v*Q.Beta)
  yh <- mu + X %*% Betah
  rss <- sum((y - yh)^2)
  R2 <- 1 - rss/s.t
  
  
  Out <- list(
    y.hat = yh,
    beta.hat = Betah,
    R2 = R2, 
    alpha = alpha
  )
  
  return(Out)
}
