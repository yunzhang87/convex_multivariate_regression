# source("first_order_opt.R")


##### QUIC path with warm start and diagonals not penalized 
### entries lambdaOmega from large to small ###
QUIC_path=function(S, lamOmega, tol = 1e-4, maxIter=1000) 
{ 
  q = dim(S)[1]; num.lam = length(lamOmega)
  out = matrix(0, q, q*num.lam)
  Sigma_prev = diag(diag(S))
  Omega_prev = 1 / Sigma_prev
  for (i in 1:num.lam) {
    lam = lamOmega[i] * (1 - diag(q))
    cov.out=QUIC(S=S, rho=lam, tol=tol, msg=0, maxIter=maxIter, X.init=Omega_prev, W.init=Sigma_prev)
    if(!is.matrix(cov.out$X)) { 
      cov.out$X=matrix(cov.out$X, nrow=q, ncol=q)
    }
      if(!is.matrix(cov.out$W)) {
        cov.out$W=matrix(cov.out$W, nrow=q, ncol=q)
      }
      Omega_prev = cov.out$X
      Sigma_prev = cov.out$W
      out[,((i-1)*q+1):(i*q)] = Omega_prev
  }
  out
}


initializer_nonconvex = function(Y, X, lamB, lamOmega, tau = 1e-2)
{
  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  sol_path = matrix(0,p+q,lamB_grid*lamOmega_grid*q)
  for (i in seq_along(lamB)){
    lam.b = lamB[i]
    C.est.ini = matrix(0,p,q)
    for (j in 1:q)
    {
      # out.lm = glmnet(X, Y[,j], family="gaussian", lambda = lamB[i], standardize = F, intercept = F)
      # C.est.ini[,j] = as.vector(out.lm$beta)
      C.est.ini[,j] = TLP(X, Y[,j], lamB[i], intercept = F)
    }
    Sigma.hat = (1/n) * crossprod(Y - X %*% C.est.ini)
    out = glasso_path_new_Rwrapper(Sigma.hat, lamOmega, tau, method=1)
    B.est.ini = C.est.ini %*% out$Omega_dc
    sol_path[1:p, ((i-1)*lamOmega_grid*q+1):(i*lamOmega_grid*q)] = B.est.ini
    sol_path[(p+1):(p+q),((i-1)*lamOmega_grid*q+1):(i*lamOmega_grid*q)] = out$Omega_dc
  }
  output = matrix(0,q,lamB_grid*lamOmega_grid*(p+q))
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid){
      output[,((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + 1):((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + p)] = t(sol_path[1:p,((i-1)*lamOmega_grid*q + (j-1)*q + 1):(((i-1)*lamOmega_grid*q) + j*q)])
      output[,((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + p + 1):((i-1)*lamOmega_grid*(p+q) + j*(p+q))] = sol_path[(p+1):(p+q),((i-1)*lamOmega_grid*q + (j-1)*q + 1):(((i-1)*lamOmega_grid*q) + j*q)]
    }
  }
  output
}

initializer_convex = function(Y, X, lamB, lamOmega)
{
  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  # Y = Y / sqrt(n); X = X / sqrt(n)
  sol_path = matrix(0,p+q,lamB_grid*lamOmega_grid*q)
  for (i in seq_along(lamB)){
    lam.b = lamB[i]
    C.est.ini = matrix(0,p,q)
    for (j in 1:q)
    {
      out.lm = glmnet(X, Y[,j], family="gaussian", lambda = lamB[i], standardize = F, intercept = F)
      C.est.ini[,j] = as.vector(out.lm$beta)
    }
    Sigma.hat = (1/n) * crossprod(Y - X %*% C.est.ini)
    # out = glasso_path_new_Rwrapper(Sigma.hat, lamOmega, tau, method=1)
    out = QUIC_path(Sigma.hat, lamOmega)
    B.est.ini = C.est.ini %*% out
    sol_path[1:p, ((i-1)*lamOmega_grid*q+1):(i*lamOmega_grid*q)] = B.est.ini
    sol_path[(p+1):(p+q),((i-1)*lamOmega_grid*q+1):(i*lamOmega_grid*q)] = out
  }
  output = matrix(0,q,lamB_grid*lamOmega_grid*(p+q))
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid){
      output[,((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + 1):((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + p)] = t(sol_path[1:p,((i-1)*lamOmega_grid*q + (j-1)*q + 1):(((i-1)*lamOmega_grid*q) + j*q)])
      output[,((i-1)*lamOmega_grid*(p+q) + (j-1)*(p+q) + p + 1):((i-1)*lamOmega_grid*(p+q) + j*(p+q))] = sol_path[(p+1):(p+q),((i-1)*lamOmega_grid*q + (j-1)*q + 1):(((i-1)*lamOmega_grid*q) + j*q)]
    }
  }
  output
}

# para = [Omega; B] \in R^{q*(p+q)}, which stacks rows of Omega and B into a vector #

grad_g_full = function(para, X, Y, XtX, YtX, YtY){
  norm_vec <- function(x) sum(x^2)
  n = dim(X)[1]; p = dim(XtX)[1]; q = dim(YtY)[1];
  ### set para to be theta and rows of B ###
  Omega = matrix(para[1:(q*q)], q, q)
  B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
  
  ### TO DO: modify gradient calculations depending on whether p,q > n or not ##
  diag.omega = diag(Omega)
  grad = rep(0,(p+q)*q)
  tmp = YtY %*% Omega - YtX %*% B
  tmp = sweep(tmp, 2, diag.omega, "/")
  grad.omega = .5 * (tmp + t(tmp))
  tmp3 = .5 * (apply(Y%*%Omega - X%*%B, 2, norm_vec) / (diag.omega^2))
  grad.omega = grad.omega - diag(tmp3)
  grad[1:(q*q)] = as.vector(grad.omega)
  tmpB = crossprod(B, XtX) - crossprod(Omega, YtX)
  tmpB = sweep(tmpB, 1, diag.omega, "/")
  grad[(q*q+1):((p+q)*q)] = as.vector(tmpB)
  # cat(para[1:9], "\n")
  # cat(grad[1:3], "\n")
  grad
}

obj_g_full = function(para, X, Y){
  norm_vec <- function(x) sum(x^2)
  q = dim(Y)[2]
  p = dim(X)[2]
  n = dim(Y)[1]
  Omega = matrix(para[1:(q*q)], q, q)
  diag.omega = diag(Omega)
  B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
  tmp = Y%*%Omega - X%*%B
  tmp.vec = apply(tmp, 2, norm_vec)
  obj_g = .5*sum(tmp.vec/diag.omega)
  # cat("omega is ", diag.omega, "\n")
  # cat("obj_g is ", obj_g, "\n")
  obj_g
}

obj_full = function(para, X, Y, lamB, lamOmega){
  norm_vec <- function(x) sqrt(sum(x^2))
  q = dim(Y)[2]
  p = dim(X)[2]
  B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
  diag.omega = diag(matrix(para[1:(q*q)], q, q))
  (obj_g_full(para, X, Y) + lamB * sum(apply(B,1,norm_vec) + lamOmega * sum(abs(para[1:(q*q)]))) - .5 * sum(log(diag.omega)))
}

obj_full_nonconvex = function(para, X, Y, lamB, lamOmega, tau) {
  norm_vec <- function(x) sqrt(sum(x^2))
  tlp <- function(x) min(abs(x)/tau, 1)
  q = dim(Y)[2]
  p = dim(X)[2]
  B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
  omega = matrix(para[1:(q*q)],q,q); tmp.diag = diag(omega); diag(omega) = .0; tmp.omega = as.vector(omega)
  penalty = lamB * sum(sapply(apply(B,1,norm_vec),tlp)) + lamOmega * sum(sapply(tmp.omega,tlp)) - .5 * sum(log(tmp.diag))
  (obj_g_full(para, X, Y) + penalty)
}

### return argmin h(x) + \|x - y\|_2^2 / (2t), ###
### where h(Omega, B) = lamB * \|B\|_{1,2} + lamOmega * \|Omega\|_{off, 1} - .5 log(diag(Omega)) ###
prox_h_full = function(y, t, lamB, lamOmega, p, q) {
  norm_vec <- function(x) sqrt(sum(x^2))
  tmp_fn <- function(x) if (x > (lamB*t)) { 1 - lamB*t/x } else { 0 }
  tmp_fn_omega <- function(x) if (abs(x) > (lamOmega*t)) { x - sign(x) * lamOmega * t } else { 0 }
  max_eps <- function(x,eps=1e-5) max(x,eps)
  B = t(matrix(y[(q*q+1):((p+q)*q)], q, p))
  row_norm <- apply(B,1,norm_vec)
  row_norm <- sapply(row_norm,tmp_fn)
  B = sweep(B,1,row_norm,"*")
  
  omega.diag = diag(matrix(y[1:(q*q)],q,q))
  omega.diag = .5 * (omega.diag + sqrt(omega.diag^2 + 2*t))
  omega.tmp = matrix(sapply(y[1:(q*q)],tmp_fn_omega),q,q)
  # diag(omega.tmp) = omega.diag
  diag(omega.tmp) = sapply(omega.diag, max_eps)
  # cat("new diag: ", diag(omega.tmp), "\n")
  c(as.vector(omega.tmp), as.vector(t(B)))
}

prox_h_full_nonconvex = function(y,t,lamB,lamOmega,tau,p,q) {
  norm_vec <- function(x) sqrt(sum(x^2))
  max_eps <- function(x,eps=1e-5) max(x,eps)
  if (lamB > 2*tau^2) {
    tmp_fn <- function(x) if (x > sqrt(2*lamB*t)) { x } else { 0 }
  } else {
    tmp_fn <- function(x) {
      if (x > (lamB*t)/(2*tau) + tau) {
        x
      } else if (x > lamB*t/tau) { 
        x - lamB*t/tau
      } else {
        0
      }
    }
  }
  if (lamOmega > 2*tau^2) {
    tmp_fn_omega <- function(x) if (abs(x) > sqrt(2*lamOmega*t)) { x } else { 0 }
  } else {
    tmp_fn_omega <- function(x) {
      if (abs(x) > (lamOmega*t)/(2*tau) + tau) {
        x
      } else if (x > lamOmega*t/tau) {
        sign(x)*(abs(x) - lamOmega*t/tau)
      } else {
        0
      }
    }
  }

  B = t(matrix(y[(q*q+1):((p+q)*q)], q, p))
  row_norm <- apply(B,1,norm_vec)
  row_norm_thred <- sapply(row_norm,tmp_fn)
  B = sweep(B,1,row_norm_thred/row_norm,"*")
  
  omega.diag = diag(matrix(y[1:(q*q)],q,q))
  # cat("old diag: ", omega.diag, "\n")
  omega.diag = .5 * (omega.diag + sqrt(omega.diag^2 + 2*t))
  omega.tmp = matrix(sapply(y[1:(q*q)],tmp_fn_omega),q,q)
  diag(omega.tmp) = sapply(omega.diag, max_eps) ## adding constraint x > eps on diag of Omega ##
  # cat("new diag", diag(omega.tmp),"\n")
  c(as.vector(omega.tmp), as.vector(t(B)))
}


# pseudo_opt_path <- function(X, Y, lamB, lamOmega, B0=NULL, Omega0 = NULL, max.iter=1e3, fast=TRUE, eps=1e-4)
# {
#   n = dim(X)[1]; p = dim(X)[2]; q = dim(Y)[2];
#   X = X /sqrt(n); Y = Y/sqrt(n);
#   YtX = crossprod(Y,X)
#   YtY = crossprod(Y); XtX = crossprod(X)
#   para = para_init = c(as.vector(diag(q)), rep(0,p*q)) 
#   if (!is.null(B0)) para_init[(q*q+1):((p+q)*q)] = as.vector(t(B0)) 
#   if (!is.null(Omega0)) para_init[1:(q*q)] = as.vector(Omega0)
  
#   obj_g = function(para) obj_g_full(para, X, Y)
#   grad_g = function(para) grad_g_full(para, X, Y, XtX, YtX, YtY)
  
#   Omega_path = array(0,c(q,q,length(lamB),length(lamOmega)))
#   B_path = array(0,c(p,q,length(lamB),length(lamOmega)))
#   for (i in seq_along(lamB)) {
#     for (j in seq_along(lamOmega)){
#       lam.b = lamB[i]; lam.omega = lamOmega[j]
#       prox_h = function(y, t) prox_h_full(y, t, lam.b, lam.omega, p, q)
#       obj = function(para) obj_full(para, X, Y, lam.b, lam.omega)
#       if (i == 1) {
#         if (fast) para = FISTA(obj_g, obj, grad_g, prox_h, x0 = para_init, max.iter=max.iter, eps=eps)
#         else para = prox_grad(obj_g, obj, grad_g, prox_h, x0 = para_init,max.iter=max.iter,eps=eps)
#       } else {
#         if (fast) para = FISTA(obj_g, obj, grad_g, prox_h, x0 = para, max.iter=max.iter,eps=eps)
#         else para = prox_grad(obj_g, obj, grad_g, prox_h, x0 = para, max.iter=max.iter,eps=eps)
#       }
#       Omega_path[,,i,j] = matrix(para[1:(q*q)], q,q)
#       B_path[,,i,j] = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
#     }
#   }
#   return(list(Omega = Omega_path, B = B_path))
# }

# B0=NULL; Omega0 = NULL; max.iter=1e3; fast=TRUE; eps=1e-3; nonconvex=TRUE; tau=NULL; warm.start=FALSE
pseudo_opt_path <- function(X, Y, lamB, lamOmega, B0=NULL, Omega0 = NULL, max.iter=1e3, fast=TRUE, eps=1e-3, nonconvex=TRUE, tau=NULL, warm.start=FALSE)
{
  n = dim(X)[1]; p = dim(X)[2]; q = dim(Y)[2];
  X = X /sqrt(n); Y = Y/sqrt(n);
  YtX = crossprod(Y,X)
  YtY = crossprod(Y); XtX = crossprod(X)
  para = para_init = c(as.vector(diag(q)), rep(0,p*q)) 
  if (!is.null(B0)) { 
    para_init[(q*q+1):((p+q)*q)] = as.vector(t(B0)) 
  }
  if (!is.null(Omega0)) { 
    para_init[1:(q*q)] = as.vector(Omega0)
  }
  if (nonconvex && is.null(tau)) {
    tau = .01
  }
  obj_g = function(para) obj_g_full(para, X, Y)
  grad_g = function(para) grad_g_full(para, X, Y, XtX, YtX, YtY)
  
  Omega_path = array(0,c(q,q,length(lamB),length(lamOmega)))
  B_path = array(0,c(p,q,length(lamB),length(lamOmega)))
  for (i in seq_along(lamB)) {
    for (j in seq_along(lamOmega)) {
      lam.b = lamB[i]; lam.omega = lamOmega[j]
      if (nonconvex) {
        prox_h = function(y, t) prox_h_full_nonconvex(y, t, lam.b, lam.omega, tau, p, q)
        obj = function(para) obj_full_nonconvex(para, X, Y, lam.b, lam.omega, tau)
      } else {
        prox_h = function(y, t) prox_h_full(y, t, lam.b, lam.omega, p, q)
        obj = function(para) obj_full(para, X, Y, lam.b, lam.omega)
      }
      if (warm.start) {
        if (i == 1 & j == 1) {
          if (fast) { 
            if (nonconvex){
              para = FISTA_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para_init, max.iter=max.iter, eps=eps) 
            } else {
              # x0 = para_init; max.iter=max.iter; eps=eps
              para = FISTA(obj_g, obj, grad_g, prox_h, x0 = para_init, max.iter=max.iter, eps=eps) 
            }
          } else {
            if (nonconvex){
              para = prox_grad_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para_init,max.iter=max.iter,eps=eps)
            } else {
              para = prox_grad(obj_g, obj, grad_g, prox_h, x0 = para_init,max.iter=max.iter,eps=eps)
            }
          }
        } else {
          if (fast) { 
            if (nonconvex){
              para = FISTA_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para, max.iter=max.iter, eps=eps) 
            } else {
              para = FISTA(obj_g, obj, grad_g, prox_h, x0 = para, max.iter=max.iter, eps=eps) 
            }
          } else {
            if (nonconvex){
              para = prox_grad_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para,max.iter=max.iter,eps=eps)
            } else {
              para = prox_grad(obj_g, obj, grad_g, prox_h, x0 = para,max.iter=max.iter,eps=eps)
            }
          }
        }
      } else {
        if (fast) { 
          if (nonconvex){
            para = FISTA_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para_init, max.iter=max.iter, eps=eps) 
          } else {
            para = FISTA(obj_g, obj, grad_g, prox_h, x0 = para_init, max.iter=max.iter, eps=eps) 
          }
        } else {
          if (nonconvex){
            para = prox_grad_nonconvex(obj_g, obj, grad_g, prox_h, x0 = para_init,max.iter=max.iter,eps=eps)
          } else {
            para = prox_grad(obj_g, obj, grad_g, prox_h, x0 = para_init,max.iter=max.iter,eps=eps)
          }
        }
      }
      Omega_path[,,i,j] = matrix(para[1:(q*q)], q,q)
      B_path[,,i,j] = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
    }
  }
  return(list(Omega.path = Omega_path, B.path = B_path))
}

# max.iter=1e3; fast=TRUE; eps=1e-4; tau =NULL;nonconvex=TRUE;cv_fold = 5;cv_method = 1

## CV uses criterion based on the loss function, not using KL loss because Omega estimate may not be PSD ##
# max.iter=1e3; fast=TRUE; eps=1e-4; tau =NULL;nonconvex=TRUE;cv_fold = 5;cv_method = 1
pseudo.cv = function(Y, X, lamB, lamOmega, max.iter=1e3, fast=TRUE, eps=1e-4,
  tau =NULL,nonconvex=TRUE,cv_fold = 5,cv_method = 1)
{
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  chunk.size = ceiling(n/cv_fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  cv.score <- matrix(0, length(lamB), length(lamOmega))
  Bpath.average <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Omegapath.average <- array(0,c(q,q,length(lamB), length(lamOmega)))
  warm_start = TRUE
  if (nonconvex) {
      warm_start = FALSE 
  }
  for (i in 1:cv_fold)
  {
    cat("CV iter ", i, " starts.... \n")
    Y.val = Y[partition[[i]],]; Y.train = Y[-partition[[i]],]
    X.val = X[partition[[i]],]; X.train = X[-partition[[i]],]
    path = pseudo_opt_path(X.train,Y.train,lamB,lamOmega,tau=tau,nonconvex=nonconvex,fast=fast,max.iter=max.iter,eps=eps, warm.start = warm_start)
    for (j in 1:length(lamB)) {
      for (k in 1:length(lamOmega)) {
        B.tmp = path$B.path[,,j,k]; Omega.tmp = path$Omega.path[,,j,k]
        # cv.score[j,k] = cv.score[j,k] + KL.dist(path$B.path[,,j,k], path$Omega.path[,,j,k], X.val, Y.val) ## may produce error b/c logdet is NA ##
        diag.omega.inv = 1 / diag(Omega.tmp)
        cv.score[j,k] = cv.score[j,k] + (1/(2*dim(X.val)[1])) * sum(diag.omega.inv * diag(crossprod(Y.val %*% Omega.tmp - X.val %*% B.tmp))) - .5 * sum(log(diag.omega.inv))
        Bpath.average[,,j,k] = Bpath.average[,,j,k] + B.tmp
        Omegapath.average[,,j,k] = Omegapath.average[,,j,k] + Omega.tmp
      }
    }
  }
  Bpath.average = Bpath.average / cv_fold; Omegapath.average = Omegapath.average / cv_fold
  cv.score = cv.score / cv_fold
  best.idx = which(cv.score==min(cv.score),arr.ind=T)
  
  B.best = Bpath.average[,,best.idx[1,1],best.idx[1,2]]
  Omega.best = Omegapath.average[,,best.idx[1,1],best.idx[1,2]]
  ## refit using selected lambdas ##
  if (cv_method == 2) {
    lamB.best = lamB[best.idx[1,1]]; lamOmega.best = lamOmega[best.idx[1,2]]
    path.final = pseudo_opt_path(X.train,Y.train,lamB.best,lamOmega.best,B0=B.best, 
      Omega0 = Omega.best, tau=tau,nonconvex=nonconvex,fast=fast,max.iter=max.iter,eps=eps)
    B.best = path.final$B.path[,,1,1]; Omega.best = path.final$Omega.path[,,1,1]
  }
  list(best.idx=best.idx,cv.score=cv.score, B.best = B.best, Omega.best = Omega.best)
}
# max.iter=1e3; fast=TRUE; eps=1e-4;Omega.initial=NULL; B.initial=NULL;nonconvex=FALSE;cv_fold = 5



