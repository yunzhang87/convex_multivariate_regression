
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

