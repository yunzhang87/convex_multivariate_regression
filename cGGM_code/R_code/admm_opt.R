source("cGGM_twostep.R")

Gene_cov<-function(p){
  sigma <- runif(p-1,0.5,1)
  covmat0 <- diag(1,p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      temp <- exp(-sum(sigma[i:(j-1)]/2))
      covmat0[i,j] <- temp
      covmat0[j,i] <- temp
    }
  }
  return(covmat0)
}

soft.threshold <- function(vec,lam){
  # Soft threshold function
  #
  # ARGUMENTS
  #	vec	: vector that soft tresholding is applied
  #	lam	: non-negative scalar- or vector-valued soft tresholding parameter. If lam is a vector, then the length of lam should be the same as that of vec
  #
  # VALUES
  #	res	: resulting vector of the same size as vec
  #
  if ( length(lam)>1 & length(lam)!=length(vec) ) {
    cat('\n ERROR: THE SIZE OF THE SECOND ARGUMENT SHOULD BE 1 OR THE SAME AS THE SIZE OF THE FIRST ARGUMENT.\n')
    return ( 0 )
  }
  vec[abs(vec)<lam] <- 0
  idx.1 <- which(vec < -lam)
  idx.2 <- which(vec > lam)
  if (length(lam)==1){
    if ( length(idx.1)>0 ) vec[idx.1]<- vec[idx.1]+lam
    if ( length(idx.2)>0 ) vec[idx.2]<- vec[idx.2]-lam
  } else if (length(lam)>1) {
    if ( length(idx.1)>0 ) vec[idx.1]<- vec[idx.1]+lam[idx.1]
    if ( length(idx.2)>0 ) vec[idx.2]<- vec[idx.2]-lam[idx.2]
  }
  return( vec )
}

soft.threshold.off.diag <- function(mat,lam){
  tmp.diag <- diag(mat)
  diag(mat) <- 0
  mat[abs(mat)<lam] <- 0
  idx.1 <- which(mat < -lam)
  idx.2 <- which(mat > lam)
  if ( length(idx.1)>0 ) mat[idx.1] <- mat[idx.1]+lam
  if ( length(idx.2)>0 ) mat[idx.2] <- mat[idx.2]-lam
  diag(mat) <- tmp.diag
  return ( mat )
}

prox_group <- function(B, lam){
  norm_vec <- function(x) sqrt(sum(x^2))
  tmp_fn <- function(x) if (x > lam) 1 - lam/x else 0
  row_norm <- apply(B,1,norm_vec)
  row_norm <- sapply(row_norm,tmp_fn)
  sweep(B,1,row_norm,"*")
}



######### optimization routine that solves multivariate regression with ####
######### D-trace loss function.                                   #########
opt_dtrace_loss <- function(Y,X,lamOmega,lamB,max.iter=1e3,rho=1) {
  norm_vec <- function(x) sqrt(sum(x^2))
  
  obj <- function(Omega, B){
    .5 * norm(crossprod(Y%*%Omega - X%*%B),"F") - sum(diag(Omega)) + 
      lamOmega * (norm(Omega,"O") - sum(abs(diag(Omega)))) + 
      lamB * sum(apply(B,1,norm_vec)) 
  }
  
  max.singular.value <- function(A){
    norm_vec <- function(x) sum(abs(x))
    row_norm <- as.matrix(apply(A,1,norm_vec))
    col_norm <- as.matrix(apply(A,2,norm_vec))
    tmp1 <- abs(t(A)) %*% row_norm
    tmp2 <- abs(A) %*% col_norm
    max(c(tmp1,tmp2))
  }
  
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2];
  X <- X/sqrt(n); Y <- Y / sqrt(n);
  gamma <- max.singular.value(cbind(Y,-X))
  Lambda <- Lambda_old <- matrix(0,n,q)
  Omega <- diag(q); B <- matrix(0,p,q);
  for (k in 1:max.iter){
    tmp <- 2*Lambda - Lambda_old
    Omega <- soft.threshold.off.diag(Omega-(crossprod(Y,tmp)-diag(q))/(rho*gamma),lamOmega/(rho*gamma))
    B <- prox_group(B+crossprod(X,tmp)/(rho*gamma),lamB/(rho*gamma))
    Lambda_old <- Lambda
    Lambda <-(Lambda + rho*(Y %*% Omega - X %*% B)) / (1 + rho)
  }
  cat("objective function is ", obj(Omega,B), "\n")
  return(list(Omega, t(solve(Omega,t(B)))))
}

### ADMM for path algorithm for dtrace loss ###
opt_dtrace_path_CD <- function(Y,X,lamOmega,lamB,max.iter=1e3,rho=1,warm.start=TRUE){
  norm_vec <- function(x) sqrt(sum(x^2))
  obj <- function(Omega, B){
    .5 * norm(crossprod(Y%*%Omega - X%*%B),"F") - sum(diag(Omega)) + 
      lamOmega * (norm(Omega,"O") - sum(abs(diag(Omega)))) + 
      lamB * sum(apply(B,1,norm_vec)) 
  }
  max.singular.value <- function(A){
    norm_vec <- function(x) sum(abs(x))
    row_norm <- as.matrix(apply(A,1,norm_vec))
    col_norm <- as.matrix(apply(A,2,norm_vec))
    tmp1 <- abs(t(A)) %*% row_norm
    tmp2 <- abs(A) %*% col_norm
    max(c(tmp1,tmp2))
  }

  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  Omega.path = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  B.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2];
  X <- X/sqrt(n); Y <- Y / sqrt(n);
  gamma <- max.singular.value(cbind(Y,-X))
  Lambda <- Lambda_old <- matrix(0,n,q)
  Omega <- diag(q); B <- matrix(0,p,q);
  
  for (i in 1:lamB_grid) {
    for (j in 1:lamOmega_grid) {
      # ## initalize (i,j)-th estimates ##
      # if (j > 1) { 
      #   Omega.ini = Omega.path[,,i,j-1]
      #   B.ini = B.path[,,i,j-1]
      # } else if (i > 1) {
      #   Omega.ini = Omega.path[,,i-1,j]
      #   B.ini = B.path[,,i-1,j]
      # } else {
      #   Omega.ini = diag(q)
      #   B.ini = matrix(0,p,q)
      # }
      lamOmega.j = lamOmega[j]
      lamB.i = lamB[i]
      for (k in 1:max.iter){
        tmp <- 2*Lambda - Lambda_old
        Omega <- soft.threshold.off.diag(Omega-(crossprod(Y,tmp)-diag(q))/(rho*gamma),lamOmega.j/(rho*gamma))
        B <- prox_group(B+crossprod(X,tmp)/(rho*gamma),lamB.i/(rho*gamma))
        Lambda_old <- Lambda
        Lambda <-(Lambda + rho*(Y %*% Omega - X %*% B)) / (1 + rho)
      }

      Omega.path[,,i,j] = Omega
      B.path[,,i,j] = B
    }
  }
  # cat("objective function is ", obj(Omega,B), "\n")

}



### path algorithm for dtrace loss with warm start ###
opt_dtrace_path_CD <- function(Y,X,lamOmega,lamB,max.iter=1e3,rho=1,warm.start=TRUE){
  lamB_grid = length(lamB)
  lamOmega_grid = length(lamOmega)
  n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
  Omega.path = array(0,dim=c(q,q,lamB_grid,lamOmega_grid))
  B.path = array(0,dim=c(p,q,lamB_grid,lamOmega_grid))
  tmp <- cGGM_twosteps(Y, X, lamB, lamOmega)
  for (i in 1:lamB_grid) {
    for (j in 1:lamOmega_grid) {
      Omega.ini = tmp$Omega.path[,,i,j]
      B.ini = tmp$C.path[,,i,j] %*% Omega.ini
      #### TO DO: implement altnerating minimization over Omega and B with warmstart ####
    }
  }
}


######### optimization routine that solves multivariate regression with ####
######### pseudo likelihood loss function. (seems to be wrong)     #########
opt_mreg_pseudo <- function(Y,X,lam,eps=.0,max.iter=1e3,rho=1,eps_abs=1e-4,eps_rel=1e-5,mu=10,eta=2){
  norm_vec <- function(x) sqrt(sum(x^2))
  
  obj <- function(Theta,B,lambda){
    obj.val = .5 * norm(crossprod(Y%*%Theta - X%*%B),"F") - sum(log(diag(Theta))) + 
      .5 * eps * norm(Theta,"F")^2  + 
      lambda * sum(apply(B,1,norm_vec)) 
    obj.val
  }
  
  max.singular.value <- function(A){
    norm_vec <- function(x) sum(abs(x))
    row_norm <- as.matrix(apply(A,1,norm_vec))
    col_norm <- as.matrix(apply(A,2,norm_vec))
    tmp1 <- abs(t(A)) %*% row_norm
    tmp2 <- abs(A) %*% col_norm
    max(c(tmp1,tmp2))
  }
  
  max.singular.value.precise <- function(A){
    eigs(crossprod(A), 1)$values
  }
  
  make_sym <- function(A){
    (A + t(A)) / 2.0
  }
  
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Y)[2];
  X <- X / sqrt(n); Y <- Y / sqrt(n);
  #gamma <- max.singular.value(cbind(Y,-X))
  gamma <- max.singular.value.precise(cbind(Y,-X))
  Lambda = Lambda_old = Lambda_delta = matrix(0,n,q)
  Theta = Theta_old = diag(q)
  B_old <- B <- matrix(0,p,q)
  B_path = list()
  Omega_path = list()
  
  ######## solution path ##########
  for (i in seq_along(lam)){
    lambda = lam[i]
    ####### admm iterations ########
    for (k in 1:max.iter){
      ### temp variables ###
      tmp <- 2*Lambda - Lambda_old
      
      ### Theta update ###
      Theta_old = Theta
      #tmp_eigen = eigen(rho * gamma * Theta - crossprod(Y, tmp))
      tmp_eigen = eigen(rho * gamma * Theta - make_sym(crossprod(Y, tmp)))
      lamlam = (tmp_eigen$values + sqrt((tmp_eigen$values)^2 + 4*(rho*gamma+eps))) / (2*(rho*gamma+eps))
      Theta = make_sym(tcrossprod(tmp_eigen$vectors, sweep(tmp_eigen$vectors, 2, lamlam,"*")))
      #Theta = tmp_eigen$vectors %*% diag(lamlam) %*% t(tmp_eigen$vectors) 
      
      ### B update ###
      B_old = B
      B <- prox_group(B + crossprod(X,tmp)/(rho*gamma),lambda/(rho*gamma))
      
      ### dual update ###
      Lambda_delta = Lambda - Lambda_old
      Lambda_old <- Lambda
      Lambda <-(Lambda + rho*(Y %*% Theta - X %*% B)) / (1 + rho)
      
      ######### stopping criterions ##########
      primal_res = norm(Lambda - Lambda_old, "F") / rho
      tmp = rho*(Y %*% (Theta-Theta_old) - X %*% (B - B_old)) + Lambda_old - Lambda + Lambda_delta
      dual_res = rho * sqrt(norm(crossprod(Y, tmp),"F")^2 + norm(crossprod(X, tmp),"F")^2)
      
      #cat("#",k,"primal residual is:",primal_res,"and dual residual is",dual_res,"\n")
      
      tmp = Y %*% Theta - X %*% B
      eps_primal = sqrt(n*q)*eps_abs + eps_rel * max(norm(tmp,"F"),norm(tmp-(Lambda-Lambda_old)/rho,"F"))
      eps_dual = sqrt((p+q)*q)*eps_abs + eps_rel * sqrt(norm(crossprod(Y,Lambda-Lambda_old),"F")^2+ norm(crossprod(X,Lambda-Lambda_old),"F")^2)
      
      if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
        cat("ADMM converges at", k, "iteration \n")
        break;
      }
      
      if ((k %% 20 == 0) && (k >= 40)) {
        #cat("objective function is ", obj(Theta,B,lambda), "\n")
        if ((primal_res / eps_primal) > mu * (dual_res / eps_dual))
          rho = rho * eta
        if (mu*(primal_res / eps_primal) <  (dual_res / eps_dual))
          rho = rho / eta
      }
    }
    
    Omega_path[[i]] = crossprod(Theta)
    B_path[[i]] = t(solve(Theta,t(B)))
    cat("objective function is ", obj(Theta,B,lambda), "\n")
  }
  return(list(Omega = Omega_path, B = B_path))
}



# ######## test #########
# require(MASS)
# n <- 200
# p <- 50
# q <- 4
# set.seed(2000)
# Sigma <- Gene_cov(q)
# Omega <- solve(Sigma)
# E <- mvrnorm(n,rep(0,q),Sigma)
# B <- matrix(0,p,q)
# B[1:3,] <- matrix(runif(3*q),3,q)
# X <- matrix(rnorm(n*p),n,p)
# Y <- X %*% B + E
# # result <- opt_dtrace_loss(Y,X,.5,.18,rho=3)
# result <- opt_mreg_pseudo(Y,X,1e-1,rho=3)
# result$Omega
# Omega
# result$B[1:10,]
# B[1:3,]
