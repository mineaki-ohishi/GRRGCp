#2020/01/19

#A using example of "f_GRRGCp.R"

n <- 50  #sample size
k <- 20  #the number of explanatory variables

X <- matrix(runif(n*k, -1, 1), n, k)  #a matrix of explanatory variables
y <- apply(X, 1, sum) + rnorm(n)      #a vector of respense variable

#prepare for obtaining an estimator
res1 <- GRRGCp(y, X, centralize=F, X.svd=F)

#get the estimator under given alpha
res2 <- get.est(res1, alpha="optimal") #optimize alpha and use the optimal alpha
#res2 <- get.est(res1, alpha=1)         #correspond to the ordinary Cp

y.hat <- res2$y.hat        #predictive value
beta.hat <- res2$beta.hat  #the estimator of regression coefficients
r2 <- res2$R2              #R-square
alpha <- res2$alpha        #the value of alpha