#!/usr/local/R/3.3.1/bin/Rscript
rm(list=ls())
sink(stdout())
args <- commandArgs(TRUE)
iter <- as.numeric(args[[1]])
print(iter)

dyn.load("../cGGM_code/lib/opt.so")
dyn.load("../cGGM_code/lib/SFGen.so")
dyn.load("../GGM_code/lib/glasso.so")

library(MASS)
library(RSpectra)
library(Matrix)
library(glmnet)
library(QUIC)
# library(MRCE)

source("../GGM_code/util_funcs.R")
source("../GGM_code/group_glasso.R")

source("TLP.R")

source("../cGGM_code/R_code/cGGM.R")
source("../cGGM_code/R_code/utils.R")
source("../cGGM_code/R_code/graph.generator.R")
source("../cGGM_code/R_code/first_order_opt.R")
source("distance.R")

## Setting
for (p in c(10)){
	for (q in c(5)){
		for (gtype in c("band", "hub", "random")){
		cat("p is ", p, "q is ", q, "graph is ", gtype, "\n")
		n <- 200
		a <- 3 ## number of nonzero rows in C
		# gtype <- "band"

		## global seed for generating all parameters ##
		set.seed(1234)
		graph <- graph.generator(q, graph = gtype)
		print(graph)
		error <- omega.generator(graph$theta, n)
		C <- matrix(0, p, q)
		C[1:a, ] <- runif(a * q, 1, 2) * sample(c(-1, 1), a * q, replace = TRUE)

		## set seed for each replicate
		set.seed(12345 + iter)

		## Data Generation
		E <- data.generator(error$omega, n)$data
		X <- matrix(rnorm(n * p), n, p)
		mu <- X %*% C
		Y <- mu + E

		## Data Analysis starts ...

		## D-trace method
		source("../cGGM_code/R_code/Dtrace_opt.R")
		lamB <- (log(p)/n) * exp(seq(from = 5, to = -1, length.out = 10))
		lamOmega <- (log(q)/n) * exp(seq(from = 4, to = 1, length.out = 10))
		dtrace.res <- dtrace.cv(Y, X, lamB, lamOmega, nonconvex = TRUE)
		print(dtrace.res$best.idx)
		C.dtrace <- dtrace.res$B.best %*% solve(dtrace.res$Omega.best)
		O.dtrace <- dtrace.res$Omega.best

		## our method based on pseudo likelihood
		source("../cGGM_code/R_code/pseudo_opt.R")
		lamB <- (log(p)/n) * exp(seq(from = 5, to = -1, length.out = 10))
		lamOmega <- (log(q)/n) * exp(seq(from = 1, to = -2, length.out = 10))
		pseudo.res <- pseudo.cv(Y, X, lamB, lamOmega)
		print(pseudo.res$best.idx)
		C.pseudo <- pseudo.res$B.best %*% solve(pseudo.res$Omega.best)
		O.pseudo <- pseudo.res$Omega.best

		## our method based on negative loglikelihood loss
		lamB = (log(p)/n) * exp(seq(from=7,to=3,length.out=10))
		lamOmega <- (log(q)/n) * exp(seq(from = 5, to = 2, length.out = 10))
		cGGM.res <- cGGM.cv_nonconvex(Y, X, lamB, lamOmega, max_newton_iter=100)
		print(cGGM.res$best.idx.nc)
		C.cGGM <- cGGM.res$B.best.nc %*% solve(cGGM.res$Omega.best.nc)
		O.cGGM <- cGGM.res$Omega.best.nc

		## MRCE group method
		source("../cGGM_code/R_code/MRCE_group/rblasso_alt.R")
		source("../cGGM_code/R_code/MRCE_group/compute.mrce.alt.R")
		source("../cGGM_code/R_code/MRCE_group/mrce.cv.R")
		source("../cGGM_code/R_code/MRCE_group/mrce.R")
		lamB <- (log(p)/n) * exp(seq(from = 3, to = -1, length.out = 10))
		lamOmega <- (log(q)/n) * exp(seq(from = 6, to = 3, length.out = 10))
		mrce.res <- mrce(X, Y, lam1.vec = lamOmega, lam2.vec = lamB, method = "cv")
		C.mrce <- mrce.res$Bhat
		O.mrce <- mrce.res$omega

		## Save results
		save.image(file = paste("result/simcv1_p", p, "_q_", q, "_graph_", gtype, "_iter", iter, ".rdata", sep = ""))
}
}
}
