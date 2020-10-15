# library(glmnet)

TLP = function(X, y, lam, family="gaussian", intercept = T, max.dc.iter = 1e2){
  small.const = 1e-8; large.const = 1 / small.const
  n = dim(X)[1] 
  y = y / sqrt(n); X = X / sqrt(n)
  X.backup = X
  nonzero_comp = integer(0)
  for (i in 1:max.dc.iter){
    nonzero_comp_old = nonzero_comp
    X[,nonzero_comp] = X.backup[,nonzero_comp] * large.const
    out = glmnet(X, y, family=family, lambda = lam, standardize = F, intercept = intercept)
    # beta.i = c(out$a0, as.vector(out$beta))
    beta.i = as.vector(out$beta)
    beta.i[nonzero_comp] = beta.i[nonzero_comp] * large.const
    nonzero_comp = which(abs(beta.i) > small.const)
    if (setequal(nonzero_comp_old, nonzero_comp)) {
      break
    }
  }
  beta.i
}


##### test #######
# set.seed(1983)
# n = 200
# p = 20
# beta =c(1,2,rep(0,p-2))
# X = matrix(rnorm(n*p), n, p)
# y = X %*% beta + rnorm(n)
# family="gaussian"; intercept = T; max.dc.iter = 1e2
# lam = (log(p) / n) * .2
# out = TLP(X, y, lam)
# out
