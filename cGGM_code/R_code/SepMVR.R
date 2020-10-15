## functions solves 1/(2n) \|Y-XC\|_F^2 + lamC * \|C\|_{1,2} 
## and its TLP version ##

# source("first_order_opt.R")
# para = [C] \in R^{q*p}, which stacks rows C into a vector #
grad_g_full = function(para, X, Y, XtX, YtX, YtY){
	n = dim(X)[1]; p = dim(XtX)[1]; q = dim(YtY)[1];
  ### set para to be rows of C ###
	C = t(matrix(para, q, p))
	as.vector(crossprod(C, XtX) - YtX)
}

obj_g_full = function(para, X, Y){
	norm_vec <- function(x) sum(x^2)
	q = dim(Y)[2]
	p = dim(X)[2]
	n = dim(Y)[1]
	C = t(matrix(para, q, p))
	obj_g = .5*norm(Y - X%*%C,"F")^2
	obj_g
}

obj_full = function(para, X, Y, lamC){
	norm_vec <- function(x) sqrt(sum(x^2))
	q = dim(Y)[2]
	p = dim(X)[2]
	C = t(matrix(para, q, p))
	penalty = lamC * sum(apply(C,1,norm_vec))
	(obj_g_full(para, X, Y) + penalty)
}

obj_full_nonconvex = function(para, X, Y, lamC, tau) {
	norm_vec <- function(x) sqrt(sum(x^2))
	tlp <- function(x) min(abs(x)/tau, 1)
	q = dim(Y)[2]
	p = dim(X)[2]
	C = t(matrix(para, q, p))
	penalty = lamC * sum(sapply(apply(C,1,norm_vec),tlp))
	(obj_g_full(para, X, Y) + penalty)
}


### return argmin h(x) + \|x - y\|_2^2 / (2t), ###
### where h(C) = lamC * \|C\|_{1,2} ###
prox_h_full = function(y, t, lamC, p, q) {
	norm_vec <- function(x) sqrt(sum(x^2))
	tmp_fn <- function(x) if (x > (lamC*t)) { 1 - lamC*t/x } else { 0 }
	C = t(matrix(y, q, p))
	row_norm <- apply(C,1,norm_vec)
	row_norm <- sapply(row_norm,tmp_fn)
	C = sweep(C,1,row_norm,"*")
	c(as.vector(t(C)))
}

prox_h_full_nonconvex = function(y,t,lamC,tau,p,q) {
	norm_vec <- function(x) sqrt(sum(x^2))
	if (lamC > 2*tau^2) {
		tmp_fn <- function(x) if (x > sqrt(2*lamC*t)) { x } else { 0 }
	} else {
		tmp_fn <- function(x) {
			if (x > (lamC*t)/(2*tau) + tau) {
				x
			} else if (x > lamC*t/tau) { 
				x - lamC*t/tau
			} else {
				0
			}
		}
	}
	C = t(matrix(y, q, p))
	row_norm <- apply(C,1,norm_vec)
	row_norm_thred <- sapply(row_norm,tmp_fn)
	C = sweep(C,1,row_norm_thred/row_norm,"*")
	c(as.vector(t(C)))
}

### nonconvex method does not work well with warm.start ###
SepMVR_opt_path <- function(X, Y, lamC, C0=NULL, 
	max.iter=1e3, fast=T, eps=1e-3, nonconvex=T, 
	tau=NULL, warm.start=F)
{
	n = dim(X)[1]; p = dim(X)[2]; q = dim(Y)[2];
	X = X /sqrt(n); Y = Y/sqrt(n);
	YtX = crossprod(Y,X)
	YtY = crossprod(Y); XtX = crossprod(X)
	para = para_init = rep(0,p*q) 
	if (!is.null(C0)) { 
		para_init = as.vector(t(C0)) 
	}
	if (nonconvex && is.null(tau)) {
		tau = .01
	}
	obj_g = function(para) obj_g_full(para, X, Y)
	grad_g = function(para) grad_g_full(para, X, Y, XtX, YtX, YtY)

	C_path = array(0,c(p,q,length(lamC)))
	for (i in seq_along(lamC)) {
		lam.c = lamC[i]
		if (nonconvex) {
			prox_h = function(y, t) prox_h_full_nonconvex(y, t, lam.c, tau, p, q)
			obj = function(para) obj_full_nonconvex(para, X, Y, lam.c, tau)
		} else {
			prox_h = function(y, t) prox_h_full(y, t, lam.c, p, q)
			obj = function(para) obj_full(para, X, Y, lam.c)
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
		C_path[,,i] = t(matrix(para, q, p))
		
	}
	return(list(C.path = C_path))
}


# max.iter=1e3; fast=TRUE; eps=1e-4; tau =NULL;nonconvex=TRUE;cv_fold = 5;cv_method = 1

