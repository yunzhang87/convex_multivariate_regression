require(MASS)
require(PDSCE)
require(expm)

norm_vec_inf <- function(X){ return(max(abs(X))) }
tr_mat <- function(X){  return(sum(diag(X)))  }
norm_vec <- function(x) sqrt(sum(x^2))


## sparsity of a square matrix excluding the diagonal elements ##
sparsity_level <- function(mat,thred=1e-10){
  d <- dim(mat)[1]
  diag(mat) = 1
  (length(which(abs(mat)>thred)) - d)/(d*(d-1))
}

loss1 <- function(omega,omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    loss <- rep(0,L)
    for (l in 1:L) {
        covmat0 <- solve(omega[[l]])
        temp_hat <- omega_hat[,(1+(l-1)*p):(l*p)]
        loss[l] = sum(apply(covmat0*temp_hat,1,sum)) - log(abs(det(covmat0%*%temp_hat))) - p
    }
    loss
}

loss2 <- function(omega,omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    loss <- rep(0,L)
    for (l in 1:L) {
        covmat0 <- solve(omega[[l]])
        temp <- covmat0 %*% omega_hat[,(1+(l-1)*p):(l*p)] - diag(p)
        loss[l] = sum(apply(temp*temp,1,sum))
    }
    loss
}

loss <- function(S,omega_hat){
    return (sum(apply(S * omega_hat,1,sum))-log(det(omega_hat)))
}

loss.like <- function(omega.hat,sigma)
{
  return ( -log(det(omega.hat)) + sum(diag(sigma%*%omega.hat)) - dim(omega.hat)[1] )
}



###### compute the average FPR ######
###### FP / (FP + TN) ###############
FP <- function(omega, omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    fp <- rep(0,L)
    for (l in 1:L){
        temp <- omega[[l]]
        temp_hat <- omega_hat[,(1+(l-1)*p):(l*p)]
        b1 <- (abs(temp) < 1e-8)
        b1 <- b1 - diag(diag(b1))
        b2 <- (abs(temp_hat) > 1e-8)
        b2 <- b2 - diag(diag(b2))
        fp[l] <- sum(apply(b1*b2,2,sum)) / sum(apply(b1,2,sum))
    }
    return(mean(fp))
}


###### compute the average FNR ######
###### FN / (TP + FN) ###############
FN <- function(omega, omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    fn <- rep(0,L)
    for (l in 1:L){
        temp <- omega[[l]]
        temp_hat <- omega_hat[,(1+(l-1)*p):(l*p)]
        b1 <- (abs(temp) > 1e-8)
        b1 <- b1 - diag(diag(b1))
        b2 <- (abs(temp_hat) < 1e-8)
        b2 <- b2 - diag(diag(b2))
        de <- sum(apply(b1,2,sum))
        if (de != 0 ){
            fn[l] <- (sum(apply(b1*b2,2,sum)) / de)
        }
    }
    return(mean(fn))
}

TP_FP <- function(omega, omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    tp <- rep(0,L)
    fp <- rep(0,L)
    Num.of.P <- rep(0,L)
    for (l in 1:L){
        temp <- omega[[l]]
        temp_hat <- omega_hat[,(1+(l-1)*p):(l*p)]
        b1 <- (abs(temp) > 1e-8)
        b1 <- b1 - diag(diag(b1))
        b2 <- (abs(temp_hat) > 1e-8)
        b2 <- b2 - diag(diag(b2))
        Num.of.P[l] <- sum(apply(b2,2,sum))
        tp[l] <- sum(apply(b1*b2,2,sum))
        fp <- Num.of.P[l] - tp[l]
    }
    return(list(TP = sum(tp), FP = sum(fp), Pos = sum(Num.of.P)))
}

TN_FN <- function(omega, omega_hat){
    p <- dim(omega_hat)[1]
    L <- dim(omega_hat)[2]/p
    tn <- rep(0,L)
    fn <- rep(0,L)
    Num.of.N <- rep(0,L)
    for (l in 1:L){
        temp <- omega[[l]]
        temp_hat <- omega_hat[,(1+(l-1)*p):(l*p)]
        b1 <- (abs(temp) < 1e-8)
        b1 <- b1 - diag(diag(b1))
        b2 <- (abs(temp_hat) < 1e-8)
        b2 <- b2 - diag(diag(b2))
        Num.of.N[l] <- sum(apply(b2,2,sum))
        tn[l] <- sum(apply(b1*b2,2,sum))
        fn <- Num.of.N[l] - tn[l]
    }
    return(list(TN = sum(tn), FN = sum(fn), Neg = sum(Num.of.N)))
}

## F1-score = 2 * (recall * precision) / (recall + precision) ##
## recall = 1 - FNR and precision = TP / (TP + FP) = Total.P * (1-FNR) / (Total.P * (1-FNR) + Total.N*FPR) #########
F1.score <- function(omega, omega_hat){
  fnr = FN(omega,omega_hat)
  ## to be implemented ##
}

performance <- function(omega,omega.hat){
    performance <- list()
    performance$fp <- FP(omega,omega.hat)
    performance$fn <- FN(omega,omega.hat)
    performance$Eloss <- loss1(omega,omega.hat)
    performance$Qloss <- loss2(omega,omega.hat)
    performance
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

rmn <- function(M, Srow, Scol){
  m=dim(Srow)[1]
  n=dim(Scol)[1]
  tmp=eigen(Srow)
  Srow.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
  tmp=eigen(Scol)
  Scol.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
  Z=matrix(rnorm(m * n), m, n)
  Srow.h %*% Z %*% Scol.h + M
}

rmn_rep <- function(N, Srow, Scol,M=NULL){
  m=dim(Srow)[1]
  n=dim(Scol)[1]
  tmp=eigen(Srow)
  Srow.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
  tmp=eigen(Scol)
  Scol.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
  #Z=matrix(rnorm(m*n*N), m, n*N)
  if (is.null(M)) { M = matrix(0,m,n) }
  Z = list()
  for (i in 1:N){
    Z[[i]] <- Srow.h %*% matrix(rnorm(m*n),m,n) %*% Scol.h + M
  }
  Z
}

check <- function(data){
  data_cleaned = list()
  for (i in seq_along(data)){
    data_cleaned[[i]] = list()  
    m <- dim(data[[i]][[1]])[1]
    n <- dim(data[[i]][[1]])[2] 
    for (j in seq_along(data[[i]])){
      cat("i,j:", dim(data[[i]][[j]])[1], dim(data[[i]][[j]])[2],"\n")
      if (dim(data[[i]][[j]])[1] == m && dim(data[[i]][[j]])[2] == n)
        { data2[[i]] = append(data2[[i]], data[[i]][[j]]) } 
      else {
        #cat("delete entry: ", i, j, "\n")
        #cat("i,j:", dim(data[[i]][[j]])[1], dim(data[[i]][[j]])[2])
      }
    }
  }
}
cor.est <- function(X){
  N <- length(X)
  m <- dim(X[[1]])[1]
  n <- dim(X[[1]])[2] 
  S.hat.row <- matrix(0,m,m)
  S.hat.col <- matrix(0,n,n)
  for (i in 1:N){
    S.hat.row <- S.hat.row + tcrossprod(X[[i]])
    S.hat.col <- S.hat.col + crossprod(X[[i]])
  }
  S.hat.row <- S.hat.row / N
  S.hat.col <- S.hat.col / N
  tmp <- diag(1/sqrt(diag(S.hat.row)))
  cor.row <- tmp %*% S.hat.row %*% tmp
  tmp <- diag(1/sqrt(diag(S.hat.col)))
  cor.col <- tmp %*% S.hat.col %*% tmp
  list(cor.row = cor.row, cor.col = cor.col,W.row = diag(S.hat.row),W.col = diag(S.hat.col))
}

###### centered correlation ###########
cor.est.new <- function(X){
  N <- length(X)
  m <- dim(X[[1]])[1]
  n <- dim(X[[1]])[2]
  X.mean <- matrix(0,m,n)
  for (i in 1:N){
    # cat(dim(X[[i]]),"\t")
    X.mean <- X.mean + data.matrix(X[[i]])
  }
  X.mean = X.mean / N
  
  S.hat.row <- matrix(0,m,m)
  S.hat.col <- matrix(0,n,n)
  for (i in 1:N){
    X[[i]] <- data.matrix(X[[i]] - X.mean)
    S.hat.row <- S.hat.row + tcrossprod(X[[i]])
    S.hat.col <- S.hat.col + crossprod(X[[i]])
  }
  S.hat.row <- S.hat.row / N
  S.hat.col <- S.hat.col / N
  tmp <- diag(1/sqrt(diag(S.hat.row)))
  cor.row <- tmp %*% S.hat.row %*% tmp
  tmp <- diag(1/sqrt(diag(S.hat.col)))
  cor.col <- tmp %*% S.hat.col %*% tmp
  list(cor.row = cor.row, cor.col = cor.col,W.row = diag(S.hat.row),W.col = diag(S.hat.col))
}


###### different time dimensions under the same group ######
cor.est.general <- function(X){
  pair_equal <- function(a, b) 
  { 
    if (a[1] == b[1] && a[2] == b[2]) {
       TRUE
    } else { FALSE }
  }
  N <- length(X)
  dim_seq = lapply(X, dim)
  dim_seq = unique(dim_seq)
  num.of.diff.time.dim = length(dim_seq)
  sample.size = rep(0,num.of.diff.time.dim)
  X.mean = list()
  for (j in 1:num.of.diff.time.dim){
    X.mean[[j]] = matrix(0,dim_seq[[j]][1], dim_seq[[j]][2])
  }
  for (i in 1:N)
  {
    for (j in 1:num.of.diff.time.dim){
      if (pair_equal(dim(X[[i]]),dim_seq[[j]])) {
        # cat("i: ",i, any(is.na(X[[i]])),"\n")
        X.mean[[j]] = X.mean[[j]] + data.matrix(X[[i]])
        sample.size[j] = sample.size[j] + 1 
      }
    }
  }
  for (j in 1:num.of.diff.time.dim){
    X.mean[[j]] = X.mean[[j]] / sample.size[j]
  }
  
  m = dim(X[[1]])[1]
  S.hat.row <- matrix(0,m,m)
  # S.hat.col <- matrix(0,n,n)
  for (i in 1:N){
    for (j in 1:num.of.diff.time.dim){
      if (pair_equal(dim(X[[i]]),dim_seq[[j]])) 
      {
          X[[i]] <- data.matrix(X[[i]] - X.mean[[j]])
          S.hat.row <- S.hat.row + tcrossprod(X[[i]])
      }
    }
  }
  S.hat.row <- S.hat.row / N
  # S.hat.col <- S.hat.col / N
  tmp <- diag(1/sqrt(diag(S.hat.row)))
  cor.row <- tmp %*% S.hat.row %*% tmp
  #tmp <- diag(1/sqrt(diag(S.hat.col)))
  #cor.col <- tmp %*% S.hat.col %*% tmp
  list(cor.row = cor.row, W.row = diag(S.hat.row))
}



hard.threshold <- function(mat,lam){
    nrow <- dim(mat)[1]
    ncol <- dim(mat)[2]
    mat[abs(mat)<lam] <- 0
    matrix(mat,nrow,ncol)
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

soft.threshold <- function(vec,lam){
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

max.singular.value <- function(A){
  norm_vec <- function(x) sum(abs(x))
  row_norm <- as.matrix(apply(A,1,norm_vec))
  col_norm <- as.matrix(apply(A,2,norm_vec))
  tmp1 <- abs(t(A)) %*% row_norm
  tmp2 <- abs(A) %*% col_norm
  max(c(tmp1,tmp2))
}

##################### for real data analysis ########################
plot.group = function(omega,omega.convex,layout=NULL,layout.grid=NULL){
  gcinfo(FALSE)  
  d <- dim(omega)[1]
  K <- dim(omega)[2]/d
  par = par(mfrow = c(2, K), pty = "s", omi=c(.3,.3,.3,.3), mai = c(.3,.3,.3,.3))
  
  tmp1 <- matrix(0,d,d)
  tmp1[abs(matrix(omega[,1:d]))>1e-10] = 1
  g = graph.adjacency(tmp1, mode="undirected", diag=FALSE)
  if (is.null(layout.grid)){
    if (is.null(layout)) { layout="auto"; layout.grid = layout.auto(g) }
    if (layout=="fruchterman.reingold") layout.grid = layout.fruchterman.reingold(g)
    if (layout=="circle") layout.grid = layout.circle(g)
  }
  for (k in 1:K){ 
    tmp1 <- matrix(0,d,d)
    tmp1[abs(matrix(omega[,((k-1)*d+1):(k*d)]))>1e-10] = 1
    g = graph.adjacency(tmp1, mode="undirected", diag=FALSE)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph nonconvex")
  }
  for (k in 1:K){ 
    tmp2 <- matrix(0,d,d)
    tmp2[abs(matrix(omega.convex[,((k-1)*d+1):(k*d)]))>1e-10] = 1
    g = graph.adjacency(tmp2, mode="undirected", diag=FALSE)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph convex")
  }
  layout.grid
}

plot.group.sub <- function(omega,ROI, layout=NULL,layout.grid=NULL,thred = 1e-3)
{
  gcinfo(FALSE)  
  d <- dim(omega)[1]
  K <- dim(omega)[2]/d
  par = par(mfrow = c(1, K), pty = "s", omi=c(.3,.3,.3,.3), mai = c(.3,.3,.3,.3))
  
  tmp1 <- matrix(0,d,d)
  tmp1[abs(matrix(omega[,1:d]))>thred] = 1
  g = graph.adjacency(tmp1, mode="undirected", diag=FALSE)
  if (is.null(layout.grid)){
    if (is.null(layout)) { layout="auto"; layout.grid = layout.auto(g) }
    if (layout=="fruchterman.reingold") layout.grid = layout.fruchterman.reingold(g)
    if (layout=="circle") layout.grid = layout.circle(g)
  }
  for (k in 1:K){ 
    omegak = omega[,((k-1)*d+1):(k*d)]
    omegak = omega[-ROI,-ROI] = 0
    tmp1 <- matrix(0,d,d)
    tmp1[abs(matrix(omega[,((k-1)*d+1):(k*d)]))>thred] = 1
    g = graph.adjacency(tmp1, mode="undirected", diag=FALSE)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph nonconvex")
  }
}


# scale = function(omega)
# {
#   d <- dim(omega)[1]
#   K <- dim(omega)[2]/d
#   for (k in 1:K){ 
#     omegak = omega[,((k-1)*d+1):(k*d)]
#     omegak = diag(1/sqrt(diag(omegak))) %*% omegak %*% diag(1/sqrt(diag(omegak)))
#     diag(omegak) = 0
#     omega[,((k-1)*d+1):(k*d)] = omegak
#   }
#   omega
# }


####### whitening matrix (to do) ###########
# require(FinCovRegularization)
## n x p matrix ##
whiten = function(X)
{
  p = dim(X)[2]
  BPHI01 = band.chol.cv(X,k.vec=NULL, method = "fast",nsplits=10,n.tr=NULL,quiet=TRUE) 
  BPHI0 = BPHI01$sigma
  de = abs (min (eigen(BPHI0)$values))+0.05
  BPHI0 = (BPHI0+de*diag(p))/(1+de)
  BPHI0INV = sqrtm(solve(BPHI0))
  return (t(as.matrix(X) %*% BPHI0INV))
}

simple_whiten = function(X)
{
  p = ncol(X); n = nrow(X)
  e = eigen(diag(p)/sqrt(n) + cor(X))
  lam = e$values
  V = e$vectors
  Omega.half.hat = diag(1/sqrt(apply(X, 2, sd))) %*% (V %*% diag(1/sqrt(lam)) %*% t(V)) %*% diag(1/sqrt(apply(X, 2, sd)))
  return (sweep(X,2,apply(X,2,mean),"-") %*% Omega.half.hat)
}


############ end util functions ##############

