KL.dist.joint <- function(Chat, Ohat, C0, Omega0, SigmaX = NULL) {
  D <- Chat - C0
  # if (is.null(SigmaX)) A <- D else A <- SigmaX %*% D
  # joint <- 0.5 * (sum((D %*% Ohat) * A) + sum(Ohat * solve(Omega0)) - log(det(Ohat)))
  dist.c <- sum(crossprod(D) * Omega0)
  dist.o <- sum(Ohat * solve(Omega0)) - log(det(Ohat)) + log(det(Omega0)) - dim(Omega0)[1]
  joint <- sum(crossprod(D) * Ohat) + dist.o
  return(c(joint, dist.c, dist.o))
}

KL.dist.path <- function(path, C0, Omega0, sigmaX = NULL, cGGM = FALSE, Cpath = FALSE) {
  if(cGGM) {
    path$B.path <- path$B.path.nc
    path$Omega.path <- path$Omega.path.nc
  }
  if(Cpath) {
    path$B.path <- path$C.path
  }
  lamB_grid <- dim(path$B.path)[3]
  lamOmega_grid <- dim(path$B.path)[4]
  dist <- matrix(0,lamB_grid, lamOmega_grid)
  for (i in 1:lamB_grid)
  {
    for (j in 1:lamOmega_grid)
    {
      Omega <- path$Omega.path[,,i,j]
      if(Cpath) {
        C <- path$B.path[,,i,j]
      } else {
        C <- path$B.path[,,i,j] %*% solve(Omega)
      }
      dist[i,j] <- KL.dist.joint(C, Omega, C0, Omega0, sigmaX)[1]
    }
  }
  best.idx <- which(dist == min(dist), arr.ind = TRUE)
  return(list(dist = dist, best.idx = best.idx))
}
  
mrce.all <- function(X, Y, lamB, lamOmega) {
  p <- dim(X)[2]
  q <- dim(Y)[2]
  Omega_path <- array(0, c(q, q, length(lamB), length(lamOmega)))
  C_path <- array(0, c(p, q, length(lamB), length(lamOmega)))
  for (i in seq_along(lamB)) {
    for(j in seq_along(lamOmega)) {
      lam2 <- lamB[i]
      lam1 <- lamOmega[j]
      res <- mrce(X, Y, lam1 = lam1, lam2 = lam2)
      C_path[, , i, j] <- res$Bhat
      Omega_path[, , i, j] <- res$omega
    }
  }
  return(list(Omega.path = Omega_path, C.path = C_path))
}

