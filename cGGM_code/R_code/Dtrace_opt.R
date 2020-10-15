# source("first_order_opt.R")
# para = [Omega; B] \in R^{q*(p+q)}, which stacks rows of Omega and B into a vector #
grad_g_full = function(para, X, Y, XtX, YtX, YtY){
	n = dim(X)[1]; p = dim(XtX)[1]; q = dim(YtY)[1];
  ### set para to be theta and rows of B ###
	Omega = matrix(para[1:(q*q)], q, q)
	B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))

  	### TO DO: modify gradient calculations depending on whether p,q > n or not ##
	grad = rep(0,(p+q)*q)
	tmp = YtY %*% Omega - YtX %*% B 
	grad.omega = (tmp + t(tmp)) / 2 - diag(q)
	# diag(grad.omega) = diag(grad.omega) / 2
	grad[1:(q*q)] = as.vector(grad.omega)
	grad[(q*q+1):((p+q)*q)] = as.vector(crossprod(B, XtX) - crossprod(Omega, YtX))
	grad
}

obj_g_full = function(para, X, Y){
	norm_vec <- function(x) sum(x^2)
	q = dim(Y)[2]
	p = dim(X)[2]
	n = dim(Y)[1]
	Omega = matrix(para[1:(q*q)], q, q)
	B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
	obj_g = .5*norm(Y%*%Omega - X%*%B,"F")^2 - sum(diag(Omega))
	obj_g
}

obj_full = function(para, X, Y, lamB, lamOmega){
	norm_vec <- function(x) sqrt(sum(x^2))
	q = dim(Y)[2]
	p = dim(X)[2]
	B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
	penalty = lamB * sum(apply(B,1,norm_vec)) + lamOmega * sum(abs(para[1:(q*q)]))
	(obj_g_full(para, X, Y) + penalty)
}

obj_full_nonconvex = function(para, X, Y, lamB, lamOmega, tau) {
	norm_vec <- function(x) sqrt(sum(x^2))
	tlp <- function(x) min(abs(x)/tau, 1)
	q = dim(Y)[2]
	p = dim(X)[2]
	B = t(matrix(para[(q*q+1):((p+q)*q)], q, p))
	omega = matrix(para[1:(q*q)],q,q); diag(omega) = .0; tmp.omega = as.vector(omega)
	penalty = lamB * sum(sapply(apply(B,1,norm_vec),tlp)) + lamOmega * sum(sapply(tmp.omega,tlp))
	(obj_g_full(para, X, Y) + penalty)
}


### return argmin h(x) + \|x - y\|_2^2 / (2t), ###
### where h(Omega, B) = lamB * \|B\|_{1,2} + lamOmega * \|Omega\|_{1} ###
prox_h_full = function(y, t, lamB, lamOmega, p, q) {
	norm_vec <- function(x) sqrt(sum(x^2))
	tmp_fn <- function(x) if (x > (lamB*t)) { 1 - lamB*t/x } else { 0 }
	tmp_fn_omega <- function(x) if (abs(x) > (lamOmega*t)) { x - sign(x) * lamOmega * t } else { 0 }
	B = t(matrix(y[(q*q+1):((p+q)*q)], q, p))
	row_norm <- apply(B,1,norm_vec)
	row_norm <- sapply(row_norm,tmp_fn)
	B = sweep(B,1,row_norm,"*")
	
	omega.diag = diag(matrix(y[1:(q*q)],q,q))
  	## force it to be positive ##
	omega.diag[omega.diag < 0] = 1e-10
	omega.tmp = matrix(sapply(y[1:(q*q)],tmp_fn_omega),q,q)
	diag(omega.tmp) = omega.diag
	c(as.vector(omega.tmp), as.vector(t(B)))
}

prox_h_full_nonconvex = function(y,t,lamB,lamOmega,tau,p,q) {
	norm_vec <- function(x) sqrt(sum(x^2))
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
	## force diag of omega to be positive ##
	omega.diag[omega.diag < 0] = 1e-10
	omega.tmp = matrix(sapply(y[1:(q*q)],tmp_fn_omega),q,q)
	diag(omega.tmp) = omega.diag
	c(as.vector(omega.tmp), as.vector(t(B)))
}

### nonconvex method does not work well with warm.start ###
Dtrace_opt_path <- function(X, Y, lamB, lamOmega, B0=NULL, Omega0 = NULL, max.iter=1e3, fast=TRUE, eps=1e-3, nonconvex=TRUE, tau=NULL, warm.start=FALSE)
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
dtrace.cv = function(Y, X, lamB, lamOmega, max.iter=1e3, fast=TRUE, eps=1e-5,
	tau =NULL,nonconvex=FALSE,cv_fold = 5,cv_method = 1)
{
	n = dim(Y)[1]; p = dim(X)[2]; q = dim(Y)[2]
	chunk.size = ceiling(n/cv_fold)
	partition <- split(1:n, ceiling(1:n/chunk.size))
	cv.score <- matrix(0, length(lamB), length(lamOmega))
	Bpath.average <- array(0, c(p, q, length(lamB), length(lamOmega)))
	Omegapath.average <- array(0,c(q,q,length(lamB), length(lamOmega)))
	if (nonconvex) {
		warm.start = FALSE
	} else {
		warm.start = TRUE
	}
	for (i in 1:cv_fold)
	{
		cat("CV iter ", i, " starts.... \n")
		Y.val = Y[partition[[i]],]; Y.train = Y[-partition[[i]],]
		X.val = X[partition[[i]],]; X.train = X[-partition[[i]],]
		path = Dtrace_opt_path(X.train,Y.train,lamB,lamOmega,tau=tau,nonconvex=nonconvex,fast=fast,max.iter=max.iter,eps=eps, warm.start=warm.start)
		for (j in 1:length(lamB)) {
			for (k in 1:length(lamOmega)) {
				B.tmp = path$B.path[,,j,k]; Omega.tmp = path$Omega.path[,,j,k]
        		# cv.score[j,k] = cv.score[j,k] + KL.dist(path$B.path[,,j,k], path$Omega.path[,,j,k], X.val, Y.val)
				cv.score[j,k] = cv.score[j,k] + (1/(2*dim(X.val)[1])) * sum((Y.val %*% Omega.tmp - X.val %*% B.tmp)^2) - sum(diag(Omega.tmp))
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
		path.final = Dtrace_opt_path(X.train,Y.train,lamB.best,lamOmega.best,B0=B.best, 
			Omega0 = Omega.best, tau=tau,nonconvex=nonconvex,fast=fast,max.iter=max.iter,eps=eps)
		B.best = path.final$B.path[,,1,1]; Omega.best = path.final$Omega.path[,,1,1]
	}
	list(best.idx=best.idx, cv.score=cv.score, B.best = B.best, Omega.best = Omega.best, B.path = Bpath.average)
}
# max.iter=1e3; fast=TRUE; eps=1e-4;Omega.initial=NULL; B.initial=NULL;nonconvex=FALSE;cv_fold = 5


