glasso_path_new_Rwrapper <- function(Sigma.hat, lam, tau,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  tau_len <- length(tau)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_len),p,p*lam_len*tau_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_len),p,p*lam_len*tau_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_new",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(tau),as.integer(tau_len),as.integer(p),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega=as.double(Omega),Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_len))
}


######## calculating path with sample covariance matrix as input #########
group_path_Rwrapper <- function(Sigma.hat,lam,nu,method=2,rho=1.0,alpha=1.5,tau=NULL,eps_abs=1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
    if (method == 1) { if (is.null(tau)) tau = .01 }
    if (method == 2) { if(is.null(tau)) tau = .1 } 
    ## tau = a * lambda in the paper. ##
    ## a = tau / lambda and R = (2 tau / lambda)^{1/2} ##
    p <- dim(Sigma.hat)[1]
    K <- dim(Sigma.hat)[2] / p
    lam_len <- length(lam)
    nu_len <- length(nu)
    if (lam[1] >= lam[lam_len]) {
      Omega <- matrix(0,p,p*K)
      for (k in 1:K)
      {
        range <- (p*(k-1)+1):(p*k)
        Omega[,range] <- diag(1/diag(Sigma.hat[,range]))
      }
      Omega <- matrix(rep(Omega,lam_len*nu_len),p,p*K*lam_len*nu_len)
    } else {
      Omega <- matrix(rep(diag(p),K*lam_len*nu_len),p,p*K*lam_len*nu_len)
      for (k in 1:K) {
          range <- (p*(k-1)+1):(p*k)
          Omega[,range] <- solve(Sigma.hat[,range]+1e-8*diag(p))
      }
    }
    out <- .C("group_nonconvex_path", as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(nu),as.integer(nu_len),as.double(tau),as.integer(p),as.integer(K),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega_L1=as.double(Omega),Omega_trL1=as.double(Omega))
    list(Omega_L1=matrix(out$Omega_L1,nrow=p,ncol=p*K*lam_len*nu_len),Omega_trL1=matrix(out$Omega_trL1,nrow=p,ncol=p*K*lam_len*nu_len))
}

group_path_constraint_Rwrapper <- function(Sigma.hat,lam,nu,R=200,tau=.01, method=2,rho=1.0,alpha=1.5,eps_abs=1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
    if (method == 1) { if (is.null(tau)) tau = .01 }
    if (method == 2) { if (is.null(tau)) tau = .1 }
    ## tau = a * lambda in the paper. ##
    ## a = tau / lambda and R = (2 tau / lambda)^{1/2} ##
    p <- dim(Sigma.hat)[1]
    K <- dim(Sigma.hat)[2] / p
    lam_len <- length(lam)
    nu_len <- length(nu)
    if (lam[1] >= lam[lam_len]) {
      Omega <- matrix(0,p,p*K)
      for (k in 1:K)
      {
        range <- (p*(k-1)+1):(p*k)
        Omega[,range] <- diag(1/diag(Sigma.hat[,range]))
      }
      Omega <- matrix(rep(Omega,lam_len*nu_len),p,p*K*lam_len*nu_len)
    } else {
      Omega <- matrix(rep(diag(p),K*lam_len*nu_len),p,p*K*lam_len*nu_len)
      for (k in 1:K) {
          range <- (p*(k-1)+1):(p*k)
          Omega[,range] <- solve(Sigma.hat[,range]+1e-8*diag(p))
      }
    }
    cat("R is ", R, "\n")
    out <- .C("group_nonconvex_constrain_path", as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(nu),as.integer(nu_len),as.double(R),as.double(tau),as.integer(p),as.integer(K),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega_L1=as.double(Omega),Omega_trL1=as.double(Omega))
    list(Omega_L1=matrix(out$Omega_L1,nrow=p,ncol=p*K*lam_len*nu_len),Omega_trL1=matrix(out$Omega_trL1,nrow=p,ncol=p*K*lam_len*nu_len))
}


######## calculating path with original data matrix as input #########
group.path <- function(data, lam, nu, method=2,rho=1,alpha=1.5,tau=1e-2,eps_abs=1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1){
    out <- list()
    num.grid <- length(lam)*length(nu)
    K = length(data)
    n.loc <- dim(data[[1]][[1]])[2]
    Sigma.loc.hat <- matrix(0,n.loc,n.loc*K)
    for (k in 1:K){
        Sigma.loc.hat[,((k-1)*n.loc+1):(k*n.loc)] <- cor.est.new(data[[k]])$cor.col
    }
    out.tmp <- group_path_Rwrapper(Sigma.loc.hat,lam,nu,method=method,rho=rho,alpha=alpha,tau=tau,eps_abs=eps_abs,eps_rel=eps_rel,N.iter=N.iter,dc.iter=dc.iter)
    out$omega.path <- list()
    out$sparsity.path <- list()
    out$omega.convex.path <- list()
    out$sparsity.convex.path <- list()
    for (k in 1:K){
        out$sparsity.path[[k]] <- rep(0,num.grid)
        out$omega.path[[k]] <- list()
        out$sparsity.convex.path[[k]] <- rep(0,num.grid)
        out$omega.convex.path[[k]] <- list()
        for (i in 1:num.grid){
            range <- ((k-1)*n.loc+1):(k*n.loc) + (i-1)*K*n.loc
            tmp <- out.tmp$Omega_trL1[,range]
            out$omega.path[[k]][[i]] <- tmp
            diag(tmp) = 0
            out$sparsity.path[[k]][i] <- length(which(abs(tmp)>1e-10))
            
            tmp <- out.tmp$Omega_L1[,range]
            out$omega.convex.path[[k]][[i]] <- tmp
            diag(tmp) = 0
            out$sparsity.convex.path[[k]][i] <- length(which(abs(tmp)>1e-10))
        }
    }
    out
}

group.cv <- function(data, lam, nu, tau=.01, num.fold=5,method=2){
    K <- length(data)
    p <- dim(data[[1]][[1]])[2]
    group.cv.tmp <- function(sigma.test,omega_path,num.of.grid){
        loss <- rep(0,num.of.grid)
        for (i in 1:num.of.grid){
            for (k in 1:K){
                loss[i] <- loss[i] + loss.like(omega_path[,((i-1)*K*p+(k-1)*p+1):((i-1)*K*p+k*p)],sigma.test[,((k-1)*p+1):(k*p)])
            }
        }
        loss/K
    }
    S.hat <- function(data) {
        S.hat <- matrix(0,p,p*K)
        for (k in 1:K) { S.hat[,((k-1)*p+1):(k*p)] <- cor.est.new(data[[k]])$cor.col }
        S.hat
    }
    n = sapply(data,length)
    partition <- list()
    for (k in 1:K) partition[[k]] <- ceiling(1:n[k]/ceiling(n[k]/num.fold))
    num.of.grid <- length(lam)*length(nu)
    data.train <- list()
    data.val <- list()
    val.idx <- list()
    omega.med <- matrix(0,p*p*K,num.fold)
    idx.best <- rep(0,num.fold)
    omega.med.convex <- matrix(0,p*p*K,num.fold)
    idx.best.convex <- rep(0,num.fold)
    for (i in 1:num.fold){
        #cat("CV fold: ",i,"\n")
        for (k in 1:K) {
            tmp <- partition[[k]]
            tmp[tmp != i] = 0
            list.tmp <- split(data[[k]],tmp)
            data.train[[k]] <- list.tmp[[1]]
            data.val[[k]] <- list.tmp[[2]]
        }
        S.train <- S.hat(data.train)
        S.val <- S.hat(data.val)
        #cat("dimension: ", dim(S.train),"\n")
        omega.path <- group_path_Rwrapper(S.train,lam,nu,method=method)
        ## nonconvex method ###
        idx.best.i <- which.min(group.cv.tmp(S.val,omega.path$Omega_trL1,num.of.grid))
        idx.best[i] <- idx.best.i
        omega.med[,i] <- as.vector(omega.path$Omega_trL1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        ### convex method ####
        idx.best.i.convex <- which.min(group.cv.tmp(S.val,omega.path$Omega_L1,num.of.grid))
        idx.best.convex[i] <- idx.best.i.convex
        omega.med.convex[,i] <- as.vector(omega.path$Omega_L1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        
    }
    omega.median <- matrix(apply(omega.med,1,median),p,p*K)
    omega.mean <- matrix(apply(omega.med,1,mean),p,p*K)
    omega.median.convex <- matrix(apply(omega.med.convex,1,median),p,p*K)
    omega.mean.convex <- matrix(apply(omega.med.convex,1,mean),p,p*K)
    list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,omega.median.convex=omega.median.convex,omega.mean.convex=omega.mean.convex,idx.best.convex=idx.best.convex)
}


## adding maximum eigenvalue constraint R = \sqrt{2 tau / lambda} ##
group.cv.constrained <- function(data, lam, nu, tau=.01, R=200, num.fold=5, method=2){
    K <- length(data)
    p <- dim(data[[1]][[1]])[2]
    group.cv.tmp <- function(sigma.test,omega_path,num.of.grid){
        loss <- rep(0,num.of.grid)
        for (i in 1:num.of.grid){
            for (k in 1:K){
                loss[i] <- loss[i] + loss.like(omega_path[,((i-1)*K*p+(k-1)*p+1):((i-1)*K*p+k*p)],sigma.test[,((k-1)*p+1):(k*p)])
            }
        }
        loss/K
    }
    S.hat <- function(data) {
        S.hat <- matrix(0,p,p*K)
        for (k in 1:K) { S.hat[,((k-1)*p+1):(k*p)] <- cor.est.new(data[[k]])$cor.col }
        S.hat
    }
    n = sapply(data,length)
    partition <- list()
    for (k in 1:K) partition[[k]] <- ceiling(1:n[k]/ceiling(n[k]/num.fold))
    num.of.grid <- length(lam)*length(nu)
    data.train <- list()
    data.val <- list()
    val.idx <- list()
    omega.med <- matrix(0,p*p*K,num.fold)
    idx.best <- rep(0,num.fold)
    omega.med.convex <- matrix(0,p*p*K,num.fold)
    idx.best.convex <- rep(0,num.fold)
    for (i in 1:num.fold){
        #cat("CV fold: ",i,"\n")
        for (k in 1:K) {
            tmp <- partition[[k]]
            tmp[tmp != i] = 0
            list.tmp <- split(data[[k]],tmp)
            data.train[[k]] <- list.tmp[[1]]
            data.val[[k]] <- list.tmp[[2]]
        }
        S.train <- S.hat(data.train)
        S.val <- S.hat(data.val)
        #cat("dimension: ", dim(S.train),"\n")
        omega.path <- group_path_constraint_Rwrapper(S.train,lam,nu,R = R,tau=tau,method=method)
        ## nonconvex method ###
        idx.best.i <- which.min(group.cv.tmp(S.val,omega.path$Omega_trL1,num.of.grid))
        idx.best[i] <- idx.best.i
        omega.med[,i] <- as.vector(omega.path$Omega_trL1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        ### convex method ####
        idx.best.i.convex <- which.min(group.cv.tmp(S.val,omega.path$Omega_L1,num.of.grid))
        idx.best.convex[i] <- idx.best.i.convex
        omega.med.convex[,i] <- as.vector(omega.path$Omega_L1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        
    }
    omega.median <- matrix(apply(omega.med,1,median),p,p*K)
    omega.mean <- matrix(apply(omega.med,1,mean),p,p*K)
    omega.median.convex <- matrix(apply(omega.med.convex,1,median),p,p*K)
    omega.mean.convex <- matrix(apply(omega.med.convex,1,mean),p,p*K)
    list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,omega.median.convex=omega.median.convex,omega.mean.convex=omega.mean.convex,idx.best.convex=idx.best.convex)
}


## real data ##
group.cv.fMRI <- function(data, lam, nu, tau=.01, num.fold=5,method=2){
    K <- length(data)
    p <- dim(data[[1]][[1]])[1]
    group.cv.tmp <- function(sigma.test,omega_path,num.of.grid){
        loss <- rep(0,num.of.grid)
        for (i in 1:num.of.grid){
            for (k in 1:K){
                loss[i] <- loss[i] + loss.like(omega_path[,((i-1)*K*p+(k-1)*p+1):((i-1)*K*p+k*p)],sigma.test[,((k-1)*p+1):(k*p)])
            }
        }
        loss/K
    }
    S.hat <- function(data) {
        S.hat <- matrix(0,p,p*K)
        for (k in 1:K) { S.hat[,((k-1)*p+1):(k*p)] <- cor.est.general(data[[k]])$cor.row }
        S.hat
    }
    n = sapply(data,length)
    partition <- list()
    for (k in 1:K) partition[[k]] <- ceiling(1:n[k]/ceiling(n[k]/num.fold))
    num.of.grid <- length(lam)*length(nu)
    data.train <- list()
    data.val <- list()
    val.idx <- list()
    omega.med <- matrix(0,p*p*K,num.fold)
    idx.best <- rep(0,num.fold)
    omega.med.convex <- matrix(0,p*p*K,num.fold)
    idx.best.convex <- rep(0,num.fold)
    for (i in 1:num.fold){
        #cat("CV fold: ",i,"\n")
        for (k in 1:K) {
            tmp <- partition[[k]]
            tmp[tmp != i] = 0
            list.tmp <- split(data[[k]],tmp)
            data.train[[k]] <- list.tmp[[1]]
            data.val[[k]] <- list.tmp[[2]]
        }
        S.train <- S.hat(data.train)
        S.val <- S.hat(data.val)
        #cat("dimension: ", dim(S.train),"\n")
        omega.path <- group_path_constraint_Rwrapper(S.train,lam,nu,tau=tau,method=method)
        ## nonconvex method ###
        idx.best.i <- which.min(group.cv.tmp(S.val,omega.path$Omega_trL1,num.of.grid))
        idx.best[i] <- idx.best.i
        cat("best index is: ", idx.best.i, "\n")
        omega.med[,i] <- as.vector(omega.path$Omega_trL1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        ### convex method ####
        idx.best.i.convex <- which.min(group.cv.tmp(S.val,omega.path$Omega_L1,num.of.grid))
        idx.best.convex[i] <- idx.best.i.convex
        omega.med.convex[,i] <- as.vector(omega.path$Omega_L1[,((idx.best.i-1)*K*p+1):(idx.best.i*K*p)])
        
    }
    omega.median <- matrix(apply(omega.med,1,median),p,p*K)
    omega.mean <- matrix(apply(omega.med,1,mean),p,p*K)
    omega.median.convex <- matrix(apply(omega.med.convex,1,median),p,p*K)
    omega.mean.convex <- matrix(apply(omega.med.convex,1,mean),p,p*K)
    list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,omega.median.convex=omega.median.convex,omega.mean.convex=omega.mean.convex,idx.best.convex=idx.best.convex)
}

group.single.cv <- function(data, lam, tau=.01, num.fold=5,method=2){
  K <- length(data)
  p <- dim(data[[1]][[1]])[2]
  glasso.cv.tmp <- function(sigma.test,omega_path,num.of.grid){
    p <- dim(sigma.test)[1]
    loss <- rep(0,num.of.grid)
    for (i in 1:num.of.grid){
      loss[i] <- loss.like(omega_path[,((i-1)*p+1):(i*p)],sigma.test)
    }
    loss
  }
  S.hat <- function(data) {
    S.hat <- matrix(0,p,p*K)
    for (k in 1:K) { S.hat[,((k-1)*p+1):(k*p)] <- cor.est.new(data[[k]])$cor.col }
    S.hat
  }
  n = sapply(data,length)
  partition <- list()
  for (k in 1:K) partition[[k]] <- ceiling(1:n[k]/ceiling(n[k]/num.fold))
  num.of.grid <- length(lam)
  data.train <- list()
  data.val <- list()
  val.idx <- list()
  omega.med <- list()
  for (k in 1:K) { omega.med[[k]] <- matrix(0,p*p,num.fold) }
  idx.best <- matrix(0,num.fold,K)
  omega.med.convex <- omega.med
  idx.best.convex <- matrix(0,num.fold,K)
  for (i in 1:num.fold){
    #cat("CV fold: ",i,"\n")
    for (k in 1:K) {
      tmp <- partition[[k]]
      tmp[tmp != i] = 0
      list.tmp <- split(data[[k]],tmp)
      data.train[[k]] <- list.tmp[[1]]
      data.val[[k]] <- list.tmp[[2]]
    }
    S.train <- S.hat(data.train)
    S.val <- S.hat(data.val)
    #cat("dimension: ", dim(S.train),"\n")
    for (k in 1:K) {
      omega.k.path <- glasso_path_new_Rwrapper(S.train[,((k-1)*p+1):(k*p)],lam,tau=c(tau),method=method)
      ## nonconvex method ###
      idx.best.k <- which.min(glasso.cv.tmp(S.val[,((k-1)*p+1):(k*p)],omega.k.path$Omega_dc,num.of.grid))
      idx.best[i,k] <- idx.best.k
      omega.med[[k]][,i] <- as.vector(omega.k.path$Omega_dc[,((idx.best.k-1)*p+1):(idx.best.k*p)])
      
      ### convex method ####
      idx.best.k <- which.min(glasso.cv.tmp(S.val[,((k-1)*p+1):(k*p)],omega.k.path$Omega,num.of.grid))
      idx.best.convex[i,k] <- idx.best.k
      omega.med.convex[[k]][,i] <- as.vector(omega.k.path$Omega[,((idx.best.k-1)*p+1):(idx.best.k*p)])
    }
  }
  omega.median <- omega.mean <- matrix(0,p,p*K)
  omega.median.convex <- omega.mean.convex <- matrix(0,p,p*K)
  for (k in 1:K){
    omega.median[,((k-1)*p+1):(k*p)] <- matrix(apply(omega.med[[k]],1,median),p,p)
    omega.mean[,((k-1)*p+1):(k*p)] <- matrix(apply(omega.med[[k]],1,mean),p,p)
    omega.median.convex[,((k-1)*p+1):(k*p)] <- matrix(apply(omega.med.convex[[k]],1,median),p,p)
    omega.mean.convex[,((k-1)*p+1):(k*p)] <- matrix(apply(omega.med.convex[[k]],1,mean),p,p)
  }
  list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,omega.median.convex=omega.median.convex,omega.mean.convex=omega.mean.convex,idx.best.convex=idx.best.convex)
}
