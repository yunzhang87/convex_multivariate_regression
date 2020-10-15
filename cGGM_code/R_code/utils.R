est_error = function(B.hat.list, B, SigmaX)
{
  error = rep(0, length(B.hat.list))
  for (i in seq_along(B.hat.list)) 
  {
    delta = B.hat.list[[i]] - B
    error[i] = sum(diag(crossprod(delta, SigmaX %*% delta)))
  }
  error
}


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

KL.dist <- function(B, Omega, X, Y)
{
  n = dim(Y)[1]
  tmp = Y %*% Omega - X %*% B
  return  (sum(solve(Omega)*crossprod(tmp))/n - log(det(Omega)))
}


KL.dist.Omega = function(Omega, Omega0)
{
  .5 * (sum(Omega*solve(Omega0)) - log(det(Omega)) + log(det(Omega0)))
}

dist.C = function(C,C0)
{
  sqrt(sum((C-C0)^2))
}

error_cGGM_C = function(result, C0)
{
  lamB_grid = dim(result$B.path)[3]
  lamOmega_grid = dim(result$B.path)[4]
  C.error = matrix(0,lamB_grid, lamOmega_grid)
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      C = result$B.path[,,i,j] %*% solve(result$Omega.path[,,i,j])
      C.error[i,j] = dist.C(C, C0)
    }
  }
  C.error
}

error_twostep_C = function(result, C0)
{
  lamC_grid = dim(result$C.path)[3]
  lamOmega_grid = dim(result$C.path)[4]
  C.error = matrix(0,lamC_grid, lamOmega_grid)
  for (i in 1:lamC_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      C = result$C.path[,,i,j]
      C.error[i,j] = dist.C(C, C0)
    }
  }
  C.error
}


error_Omega = function(result, Omega0)
{
  lamB_grid = dim(result$Omega.path)[3]
  lamOmega_grid = dim(result$Omega.path)[4]
  Omega.error = matrix(0,lamB_grid, lamOmega_grid)
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      Omega = result$Omega.path[,,i,j]
      Omega.error[i,j] = KL.dist.Omega(Omega, Omega0)
    }
  }
  Omega.error
}

error_cGGM = function(result, Omega0, C0)
{
  lamB_grid = dim(result$B.path)[3]
  lamOmega_grid = dim(result$B.path)[4]
  error = matrix(0,lamB_grid, lamOmega_grid)
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      Omega = result$Omega.path[,,i,j]
      C = result$B.path[,,i,j] %*% solve(Omega)
      error[i,j] = sum(Omega*solve(Omega0)) - log(det(Omega)) + log(det(Omega0)) + sum(result$Omega.path[,,i,j] * crossprod(C - C0))
    }
  }
  error
}

error_twostep = function(result, Omega0, C0)
{
  lamB_grid = dim(result$C.path)[3]
  lamOmega_grid = dim(result$C.path)[4]
  error = matrix(0,lamB_grid, lamOmega_grid)
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      Omega = result$Omega.path[,,i,j]
      C = result$C.path[,,i,j]
      error[i,j] = sum(Omega*solve(Omega0)) - log(det(Omega)) + log(det(Omega0)) + sum(result$Omega.path[,,i,j] * crossprod(C - C0))
    }
  }
  error
}

KL.dist.joint <- function(Chat, Ohat, C0, Omega0, SigmaX = NULL) {
  D <- Chat - C0
  if (is.null(SigmaX)) A <- D else A <- SigmaX %*% D
  joint <- 0.5 * (sum((D %*% Ohat) * A) + sum(Ohat * solve(Omega0)) - log(det(Ohat)))
  dist.c <- sum(crossprod(D) * Omega0)
  dist.o <- sum(Ohat * solve(Omega0)) - log(det(Ohat))
  return(c(joint, dist.c, dist.o))
}

