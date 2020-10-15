rm(list=ls())
source("../first_order_opt.R")
source("mrce.R")
source("../utils.R")
require(MASS)
require(RSpectra)
n <- 200
p <- 100
q <- 3

set.seed(2000)
Sigma <- Gene_cov(q)
Omega <- solve(Sigma)
E <- mvrnorm(n,rep(0,q),Sigma)
C <- matrix(0,p,q)
C[1:3,] <- matrix(3,3,q)
X <- matrix(rnorm(n*p),n,p)
Y <- X %*% C + E

n = dim(Y)[1]
# lamB = (log(p)/n) * exp(seq(from=1,to=-4,length.out=40))
# lamOmega = (log(p)/n) * exp(seq(from=2,to=0,length.out=10))


lamC = .01
lamOmega = .0
out = mrce(X, Y, lam1 = lamOmega, lam2 = lamC, method="single", silent=F)
out$Bhat[1:10,]

lamC = c(.01,.009)
lamOmega = c(.0,.001)
out = mrce(X, Y, lam1.vec = lamOmega, lam2.vec = lamC, method="cv", silent=T)
