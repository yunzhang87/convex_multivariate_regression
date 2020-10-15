######### function get a path of tlp and lasso solution for #########
######### precision matrix                                  #########
library(QUIC)

TLP_matrix <- function(S, Lambda, tau=1e-2,
  tol=1e-4, maxIter=1000, max.dc.iter=20,weighted=TRUE) {
  num.of.grid <- length(Lambda)
  p <- dim(S)[1]
  lambda.max <- max(abs(S))
  rho <- matrix(1,p,p)
  if (weighted) {
    diag.S <- diag(S)
    rho <- matrix(sqrt(diag.S)) %*% sqrt(diag.S)
    lambda.max <- max(abs(S)/rho)
  }
  diag(rho) <- 0
  sol.con <- QUIC(S,rho*lambda.max,path = Lambda,
    tol=tol,msg = 0,maxIter=maxIter)
  X.lasso <- sol.con$X
  W.lasso <- sol.con$W
  X.tlp <- X.lasso
  W.tlp <- W.lasso
  for (i in 1:num.of.grid) {
    lambda <- Lambda[i]*lambda.max
    for (dc.iter in 1:max.dc.iter) {
      rho.dc <- rho
      rho.dc[abs(X.tlp[,,i])>tau] <- 0
      sol <- QUIC(S,rho.dc*lambda,tol=tol,msg=0,maxIter=maxIter,
        X.init=X.tlp[,,i],W.init=W.tlp[,,i])
      tmp <- sol$X
      if (mean((tmp-X.tlp[,,i])^2) < 1e-8) { 
        break 
      } else {
        X.tlp[,,i] <- tmp
        W.tlp[,,i] <- sol$W
      }
    }
  }
  results <- list()
  results$X.lasso <- X.lasso
  results$X.tlp <- X.tlp
  return(results)
}


