initializer_path_tlp = function(Y, X, lamB, lamOmega, tau = 1e-2)
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

initializer_path_lasso = function(Y, X, lamB, lamOmega, tau =1e-2)
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

post_processing <- function(Y, X, B.hat, lamOmega, gamma=.01) {
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  C.oracle = matrix(0,p,q)
  norm.2 = function(x) sqrt(sum(x^2))
  S.hat = which(apply(B.hat, 1, norm.2) > (sqrt(q) * gamma))
  if (length(S.hat) == 0) {
    Sigma.hat = crossprod(Y) / n
  } else {
    tmp = crossprod(X[,S.hat]) + 1e-8 * diag(length(S.hat))
    # if (length(S.hat) >= n) tmp = tmp + 1e-8 * diag(length(S.hat))
    C.hat = solve(tmp, crossprod(X[,S.hat],Y))
    C.oracle[S.hat, ] = C.hat
    Sigma.hat = crossprod(Y - X %*% C.oracle) / n
  }
  Omega.hat = glasso_path_new_Rwrapper(Sigma.hat, lamOmega, gamma, method=1)$Omega_dc
  B.hat = C.oracle %*% Omega.hat 
  list(B.hat = B.hat, Omega.hat = Omega.hat)
}


cGGM_nonconvex <- function(Y, X, lamB, lamOmega, gamma = .01, eps_abs=1e-7, 
                   eps_rel=1e-8, newton_tol=1e-5, mu=10.0, tau=2.0, rho=1.0, 
                   alpha=.4, max_newton_iter=10, max_admm_iter=1e3,
                   max_dc_iter = 5,ini="lasso", post_process=T, nonconvex = 1)
{
  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  if (ini=="lasso") {
    sol_path_ini = initializer_path_lasso(Y, X, lamB, lamOmega) # matrix(0,p+q,lamB_grid*lamOmega_grid*q)
    sol_path_ini_lasso = sol_path_ini
  } else {
    sol_path_ini = initializer_path_tlp(Y, X, lamB, lamOmega)
  }
  sol_path_nc = matrix(0,p+q,lamB_grid*lamOmega_grid*q)
  Z.initial = sol_path_ini[,1:q]
  res = .C("cGGM_nonconvex", as.double(t(Y)), as.double(t(X)), as.double(t(Z.initial)), as.double(lamB), 
           as.double(lamOmega), gamma = as.double(gamma), sol_path = as.double(sol_path_ini), 
           sol_path_nc = as.double(sol_path_nc), as.integer(p), as.integer(q), 
           as.integer(n), as.integer(lamB_grid), as.integer(lamOmega_grid), as.double(eps_abs), 
           as.double(eps_rel), as.double(newton_tol), as.double(mu), as.double(tau), as.double(rho), 
           as.double(alpha), as.integer(max_newton_iter), as.integer(max_admm_iter), 
           as.integer(max_dc_iter), as.integer(nonconvex))
  sol_path = t(matrix(res$sol_path, nrow = q, lamB_grid*lamOmega_grid*(p+q)))
  sol_path_nc = t(matrix(res$sol_path_nc, nrow = q, lamB_grid*lamOmega_grid*(p+q)))
  Omega.path = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  #C.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  B.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))

  Omega.path.ini = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  B.path.ini = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))

  Omega.path.nc = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  #C.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  B.path.nc = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  for (i in 1:lamB_grid){
    for (j in 1:lamOmega_grid){
      Omega.path[,,i,j] = sol_path[(p+q)*((i-1)*lamOmega_grid+j-1) + (p+1):(p+q),]
      B.path[,,i,j] = sol_path[(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p,]

      Omega.path.ini[,,i,j] = t(sol_path_ini_lasso[,(p+q)*((i-1)*lamOmega_grid+j-1) + (p+1):(p+q)])
      B.path.ini[,,i,j] = t(sol_path_ini_lasso[,(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p])

      if (post_process) {
        B.tmp = sol_path_nc[(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p,]
        out = post_processing(Y, X, B.tmp, lamOmega[j], gamma=gamma)
        Omega.path.nc[,,i,j] = out$Omega.hat 
        B.path.nc[,,i,j] = out$B.hat 
      } else {
        Omega.path.nc[,,i,j] = sol_path_nc[(p+q)*((i-1)*lamOmega_grid+j-1) + (p+1):(p+q),]
        B.path.nc[,,i,j] = sol_path_nc[(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p,]
      }
      #C.path[,,i,j] = B.path[,,i,j] %*% solve(Omega.path[,,i,j])
    }
  }
  
    list(Omega.path = Omega.path, B.path = B.path,
       Omega.path.nc = Omega.path.nc, B.path.nc = B.path.nc, 
       Omega.path.ini = Omega.path.ini, B.path.ini = B.path.ini)
}

cGGM.cv_nonconvex = function(Y, X, lamB, lamOmega, gamma = .01, 
                   Omega.initial=NULL, B.initial=NULL, eps_abs=1e-7, 
                   eps_rel=1e-8, newton_tol=1e-5, mu=10.0, tau=2.0, rho=1.0, 
                   alpha=.4, max_newton_iter=10, max_admm_iter=1e3,
                   max_dc_iter = 5, cv_fold = 5, ini="lasso", post_process=T, 
                   nonconvex = 1)
{
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  chunk.size = ceiling(n/cv_fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  cv.score <- matrix(0, length(lamB), length(lamOmega))
  cv.score.ini <- matrix(0, length(lamB), length(lamOmega))
  cv.score.nc <- matrix(0, length(lamB), length(lamOmega))
  Bpath.average <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Bpath.average.ini <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Bpath.average.nc <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Omegapath.average <- array(0,c(q,q,length(lamB), length(lamOmega)))
  Omegapath.average.nc <- array(0,c(q,q,length(lamB), length(lamOmega)))
  Omegapath.average.ini <- array(0,c(q,q,length(lamB), length(lamOmega)))
  for (i in 1:cv_fold)
  {
    cat("CV iter ", i, " starts.... \n")
    Y.val = Y[partition[[i]],]; Y.train = Y[-partition[[i]],]
    X.val = X[partition[[i]],]; X.train = X[-partition[[i]],]
    path = cGGM_nonconvex(Y.train,X.train,lamB,lamOmega,gamma=gamma,
               eps_abs=eps_abs, eps_rel=eps_rel, newton_tol=newton_tol, 
               mu=mu, tau=tau, rho=rho, alpha=alpha, 
               max_newton_iter=max_newton_iter, max_admm_iter=max_admm_iter,
               max_dc_iter = max_dc_iter,ini=ini, post_process=post_process, nonconvex = nonconvex)
    for (j in 1:length(lamB)){
      for (k in 1:length(lamOmega)){
        n.val = dim(Y.val)[1]
        B.tmp = path$B.path[,,j,k]; Omega.tmp = path$Omega.path[,,j,k]
        B.tmp.nc = path$B.path.nc[,,j,k]; Omega.tmp.nc = path$Omega.path.nc[,,j,k]
        B.tmp.ini = path$B.path.ini[,,j,k]; Omega.tmp.ini = path$Omega.path.ini[,,j,k]
        cv.score[j,k] = cv.score[j,k] + sum(solve(Omega.tmp) * crossprod(Y.val %*% Omega.tmp - X.val %*% B.tmp)) - n.val * log(det(Omega.tmp))
        cv.score.nc[j,k] = cv.score.nc[j,k] + sum(solve(Omega.tmp.nc) * crossprod(Y.val %*% Omega.tmp.nc - X.val %*% B.tmp.nc)) - n.val * log(det(Omega.tmp.nc))
        cv.score.ini[j,k] = cv.score.ini[j,k] + sum(solve(Omega.tmp.ini) * crossprod(Y.val %*% Omega.tmp.ini - X.val %*% B.tmp.ini)) - n.val * log(det(Omega.tmp.ini))
        Bpath.average[,,j,k] = Bpath.average[,,j,k] + B.tmp
        Omegapath.average[,,j,k] = Omegapath.average[,,j,k] + Omega.tmp
        Bpath.average.nc[,,j,k] = Bpath.average.nc[,,j,k] + B.tmp.nc
        Omegapath.average.nc[,,j,k] = Omegapath.average.nc[,,j,k] + Omega.tmp.nc
        Bpath.average.ini[,,j,k] = Bpath.average.ini[,,j,k] + B.tmp.ini
        Omegapath.average.ini[,,j,k] = Omegapath.average.ini[,,j,k] + Omega.tmp.ini
        # cv.score[j,k] = cv.score[j,k] + KL.dist(B.tmp, Omega.tmp, X.val, Y.val)
        # cv.score.nc[j,k] = cv.score.nc[j,k] + KL.dist(path$B.path.nc[,,j,k], path$Omega.path.nc[,,j,k], X.val, Y.val)
      }
    }
  }
  cv.score = cv.score / cv_fold
  cv.score.nc = cv.score.nc / cv_fold
  cv.score.ini = cv.score.ini / cv_fold
  best.idx = which(cv.score==min(cv.score),arr.ind=T)
  best.idx.nc = which(cv.score.nc==min(cv.score.nc),arr.ind=T)
  best.idx.ini = which(cv.score.ini==min(cv.score.ini),arr.ind=T)
  
  B.best = Bpath.average[,,best.idx[1,1],best.idx[1,2]]/cv_fold
  Omega.best = Omegapath.average[,,best.idx[1,1],best.idx[1,2]]/cv_fold
  B.best.nc = Bpath.average.nc[,,best.idx.nc[1,1],best.idx.nc[1,2]]/cv_fold
  Omega.best.nc = Omegapath.average.nc[,,best.idx.nc[1,1],best.idx.nc[1,2]]/cv_fold
  B.best.ini = Bpath.average.ini[,,best.idx.ini[1,1],best.idx.ini[1,2]]/cv_fold
  Omega.best.ini = Omegapath.average.ini[,,best.idx.ini[1,1],best.idx.ini[1,2]]/cv_fold
  list(B.best = B.best, Omega.best = Omega.best, B.best.nc = B.best.nc, B.best.ini = B.best.ini, Omega.best.nc = Omega.best.nc, Omega.best.ini = Omega.best.ini, best.idx=best.idx,best.idx.nc=best.idx.nc, best.idx.ini = best.idx.ini, cv.score=cv.score,cv.score.nc=cv.score.nc, cv.score.ini=cv.score.ini)
}

# Omega.initial=NULL; B.initial=NULL; eps_abs=1e-7; 
# eps_rel=1e-8; newton_tol=1e-5; mu=10.0; tau=2.0; rho=1.0; 
# alpha=.4; max_newton_iter=1e2; max_admm_iter=5e3
cGGM <- function(Y, X, lamB, lamOmega, gamma = .01, 
                   Omega.initial=NULL, B.initial=NULL, eps_abs=1e-7, 
                   eps_rel=1e-8, newton_tol=1e-5, mu=10.0, tau=2.0, rho=1.0, 
                   alpha=.4, max_newton_iter=2e2, max_admm_iter=5e3,
                   max_dc_iter = 5, nonconvex=1)
{
  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  if (is.null(Omega.initial)) Omega.initial = diag(q)
  if (is.null(B.initial)) B.initial = matrix(0, p, q)
  Z.initial = rbind(B.initial,Omega.initial)
  sol_path = matrix(0,p+q,lamB_grid*lamOmega_grid*q)
  sol_path_nc = NULL
  if (nonconvex == 1) sol_path_nc = matrix(0,p+q,lamB_grid*lamOmega_grid*q)
  res = .C("cGGM", as.double(t(Y)), as.double(t(X)), as.double(t(Z.initial)), as.double(lamB), 
           as.double(lamOmega), gamma = as.double(gamma), sol_path = as.double(sol_path), 
           sol_path_nc = as.double(sol_path_nc), as.integer(p), as.integer(q), 
           as.integer(n), as.integer(lamB_grid), as.integer(lamOmega_grid), as.double(eps_abs), 
           as.double(eps_rel), as.double(newton_tol), as.double(mu), as.double(tau), as.double(rho), 
           as.double(alpha), as.integer(max_newton_iter), as.integer(max_admm_iter), 
           as.integer(max_dc_iter), as.integer(nonconvex))
  sol_path = t(matrix(res$sol_path, nrow = q, lamB_grid*lamOmega_grid*(p+q)))
  if (nonconvex == 1) sol_path_nc = t(matrix(res$sol_path_nc, nrow = q, lamB_grid*lamOmega_grid*(p+q)))
  Omega.path = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  #C.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  B.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  if (nonconvex == 1) Omega.path.nc = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  #C.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  if (nonconvex == 1) B.path.nc = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  for (i in 1:lamB_grid){
    for (j in 1:lamOmega_grid){
      Omega.path[,,i,j] = sol_path[(p+q)*((i-1)*lamOmega_grid+j-1) + (p+1):(p+q),]
      B.path[,,i,j] = sol_path[(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p,]
      if (nonconvex == 1) Omega.path.nc[,,i,j] = sol_path_nc[(p+q)*((i-1)*lamOmega_grid+j-1) + (p+1):(p+q),]
      if (nonconvex == 1) B.path.nc[,,i,j] = sol_path_nc[(p+q)*((i-1)*lamOmega_grid+j-1) + 1:p,]
      #C.path[,,i,j] = B.path[,,i,j] %*% solve(Omega.path[,,i,j])
    }
  }
  if (nonconvex == 1) {
    list(Omega.path = Omega.path, B.path = B.path,
       #, C.path = C.path 
       Omega.path.nc = Omega.path.nc, B.path.nc = B.path.nc)
  } else {
    list(Omega.path = Omega.path, B.path = B.path)
  }
}


cGGM.cv = function(Y, X, lamB, lamOmega, gamma = .01, 
                   Omega.initial=NULL, B.initial=NULL, eps_abs=1e-7, 
                   eps_rel=1e-8, newton_tol=1e-5, mu=10.0, tau=2.0, rho=1.0, 
                   alpha=.4, max_newton_iter=1e2, max_admm_iter=5e3,
                   max_dc_iter = 5, nonconvex=1, cv_fold = 5)
{
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  chunk.size = ceiling(n/cv_fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  cv.score <- matrix(0, length(lamB), length(lamOmega))
  cv.score.nc <- matrix(0, length(lamB), length(lamOmega))
  Bpath.average <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Omegapath.average <- array(0,c(q,q,length(lamB), length(lamOmega)))
  Bpath.average.nc <- array(0, c(p,q, length(lamB), length(lamOmega)))
  Omegapath.average.nc <- array(0,c(q,q,length(lamB), length(lamOmega)))
  for (i in 1:cv_fold)
  {
    cat("CV iter ", i, " starts.... \n")
    Y.val = Y[partition[[i]],]; Y.train = Y[-partition[[i]],]
    X.val = X[partition[[i]],]; X.train = X[-partition[[i]],]
    path = cGGM(Y.train,X.train,lamB,lamOmega,gamma=gamma,
               Omega.initial=Omega.initial, B.initial=B.initial, 
               eps_abs=eps_abs, eps_rel=eps_rel, newton_tol=newton_tol, 
               mu=mu, tau=tau, rho=rho, alpha=alpha, 
               max_newton_iter=max_newton_iter, max_admm_iter=max_admm_iter,
               max_dc_iter = max_dc_iter, nonconvex=nonconvex)
    for (j in 1:length(lamB)){
      for (k in 1:length(lamOmega)){
        n.val = dim(Y.val)[1]
        B.tmp = path$B.path[,,j,k]; Omega.tmp = path$Omega.path[,,j,k]
        if (nonconvex == 1) B.tmp.nc = path$B.path.nc[,,j,k]; Omega.tmp.nc = path$Omega.path.nc[,,j,k]
        cv.score[j,k] = cv.score[j,k] + sum(solve(Omega.tmp) * crossprod(Y.val %*% Omega.tmp - X.val %*% B.tmp)) - n.val * log(det(Omega.tmp))
        if (nonconvex == 1) cv.score.nc[j,k] = cv.score.nc[j,k] + sum(solve(Omega.tmp.nc) * crossprod(Y.val %*% Omega.tmp.nc - X.val %*% B.tmp.nc)) - n.val * log(det(Omega.tmp.nc))
        Bpath.average[,,j,k] = Bpath.average[,,j,k] + B.tmp
        Omegapath.average[,,j,k] = Omegapath.average[,,j,k] + Omega.tmp
        if (nonconvex == 1)  {
          Bpath.average.nc[,,j,k] = Bpath.average.nc[,,j,k] + B.tmp.nc
          Omegapath.average.nc[,,j,k] = Omegapath.average.nc[,,j,k] + Omega.tmp.nc
        }
        # cv.score[j,k] = cv.score[j,k] + KL.dist(B.tmp, Omega.tmp, X.val, Y.val)
        # cv.score.nc[j,k] = cv.score.nc[j,k] + KL.dist(path$B.path.nc[,,j,k], path$Omega.path.nc[,,j,k], X.val, Y.val)
      }
    }
  }
  cv.score = cv.score / cv_fold
  cv.score.nc = cv.score.nc / cv_fold
  best.idx = which(cv.score==min(cv.score),arr.ind=T)
  best.idx.nc = which(cv.score.nc==min(cv.score.nc),arr.ind=T)
  
  B.best = Bpath.average[,,best.idx[1,1],best.idx[1,2]]/cv_fold
  Omega.best = Omegapath.average[,,best.idx[1,1],best.idx[1,2]]/cv_fold
  B.best.nc = Bpath.average.nc[,,best.idx.nc[1,1],best.idx.nc[1,2]]/cv_fold
  Omega.best.nc = Omegapath.average.nc[,,best.idx.nc[1,1],best.idx.nc[1,2]]/cv_fold
  list(B.best = B.best, Omega.best = Omega.best, B.best.nc = B.best.nc, Omega.best.nc = Omega.best.nc, best.idx=best.idx,best.idx.nc=best.idx.nc,cv.score=cv.score,cv.score.nc=cv.score.nc)
}

# gamma = .01; 
# Omega.initial=NULL; B.initial=NULL; eps_abs=1e-7; 
# eps_rel=1e-8; newton_tol=1e-5; mu=10.0; tau=2.0; rho=1.0; 
# alpha=.4; max_newton_iter=1e2; max_admm_iter=5e3;
# max_dc_iter = 5; nonconvex=0; cv_fold = 5


