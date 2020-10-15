## minimize g(x) + h(x) where g(x) is smooth and h(x) is nonsmooth but proximable ##

# eta = .8; t = 1e-1; eps=1e-4
FISTA <- function(obj_g, obj, grad_g, prox_h, x0=NULL, dim=NULL, max.iter = 1e4, beta = .8, t = 1e-1, eps=1e-4){
  if (is.null(x0)) {
    if (is.null(dim)) {
      stop("Please specify either an initial solution or input dimension")
    }
    x0 = rep(0,dim)
  } else {
    if (is.null(dim)) { 
      dim = length(x0) 
    } else { 
      if (dim != length(x0)) {
        stop("dimension of initial solution must have dimension dim") 
      }
    }
  }
  x_prev = x_cur = z = x0 
  delta_x = rep(0,dim)
  for (k in 1:max.iter){
    delta_x = x_cur - x_prev
    z = x_cur + delta_x * (k-1)/(k+2)
    x_prev = x_cur
    x_cur = prox_h(z-t*grad_g(z), t)
    delta_x = x_cur - z
    # cat("step size is", t, "\n")
    # cat("objective function is ", obj(x_cur), "\n")
    while (obj_g(x_cur) > (obj_g(z) + crossprod(grad_g(z),delta_x) + sum(delta_x^2)/(2*t))){
     t = beta * t
     x_cur = prox_h(z-t*grad_g(z), t)
     delta_x = x_cur -z
   }

   ####### check convergence #######
   delta_x = x_cur - prox_h(x_cur-grad_g(x_cur), 1)
   if (max(abs(delta_x)) < eps) {
    # cat("Converges at", k, "iteration.  \n")
    # cat("objective function is ", obj(x_cur), "\n")
    break; 
  }
  if (k == max.iter) {
    # cat("max number of iteration reached! \n")
    # cat("objective function is ", obj(x_cur), "\n")
  }
}
x_cur
}

## minimize g(x) + h(x) where g(x) is smooth and h(x) is nonsmooth but proximable ##
prox_grad <- function(obj_g, obj, grad_g, prox_h, x0=NULL, dim=NULL, max.iter = 2e3, beta = .8, t = 1e-2,eps = 1e-5){
  if (is.null(x0)) {
    if (is.null(dim)) stop("Please specify either an initial solution or input dimension")
      x0 = rep(0,dim)
  }
  else {
    if (is.null(dim)) { dim = length(x0) }
    else { if (dim != length(x0)) stop("dimension of initial solution must have dimension dim") }
  }
  x_prev = x_cur = z = x0
  delta_x = rep(0,dim)
  for (k in 1:max.iter) {
    x_prev = x_cur
    x_cur = prox_h(x_prev-t*grad_g(x_prev), t)
    delta_x = x_cur - x_prev
    # cat("step size is", t, "\n")
    # cat("objective function is ", obj(x_prev), "\n")
    while (obj_g(x_cur) > (obj_g(x_prev) + crossprod(grad_g(x_prev), delta_x) + sum(delta_x^2) / (2*t))) {
      t = t * beta
      x_cur = prox_h(x_prev-t*grad_g(x_prev), t)
      delta_x = x_cur - x_prev
    }
    ####### check convergence #######
    delta_x = x_cur - prox_h(x_cur-grad_g(x_cur), 1)
    if (max(abs(delta_x)) < eps) {
     #cat("Converges at", k, "iteration.  \n")
     #cat("objective function is ", obj(x_cur), "\n")
      ##cat(x_cur[1:9], "\n")
      ##cat(grad_g(x_cur)[1:9], "\n")
      break; 
    }
    if (k == max.iter) {
     #cat("max number of iteration reached! \n")
     #cat("objective function is ", obj(x_cur), "\n")
      # cat(x_cur[1:9], "\n")
      # cat(grad_g(x_cur)[1:9], "\n")
    }
  }
  
  x_cur
}


FISTA_nonconvex <- function(obj_g, obj, grad_g, prox_h, x0=NULL, dim=NULL, max.iter = 1e4, beta = .8, t = 1e-1, eps=1e-4){
  if (is.null(x0)) {
    if (is.null(dim)) stop("Please specify either an initial solution or input dimension")
      x0 = rep(0,dim)
  }
  else {
    if (is.null(dim)) { dim = length(x0) }
    else { if (dim != length(x0)) stop("dimension of initial solution must have dimension dim") }
  }
  x_prev = x_cur = z = x0 
  delta_x = rep(0,dim)
  for (k in 1:max.iter){
    delta_x = x_cur - x_prev
    z = x_cur + delta_x * (k-1)/(k+2)
    # cat("k is: ",k, " and z is: ", z[1:9],"\n")
    x_prev = x_cur
    x_cur = prox_h(z-t*grad_g(z), t)
    delta_x = x_cur - z
    
    # cat("step size is", t, "\n")
    # cat("objective function is ", obj(x_cur), "\n")
    
    while (obj_g(x_cur) > (obj_g(z) + crossprod(grad_g(z),delta_x) + sum(delta_x^2)/(2*t))){
     t = beta * t
     x_cur = prox_h(z-t*grad_g(z), t)
     delta_x = x_cur - z
     # cat("step size is ", t, "\n")
     # cat((obj_g(z) + crossprod(grad_g(z),delta_x) + sum(delta_x^2)/(2*t)),"\n")
     # cat(x_cur, "\n")
    # cat(obj_g(x_cur), "\n")
    #  if (obj_g(x_cur) == Inf){ stop() }
   }

    ####### check convergence #######
   # cat("current iterate: ", x_cur, "\n")
   # cat("previous iterate: ", x_prev, "\n")
   # cat("current obj: ", obj(x_cur), "\n")
   # cat("previous obj: ", obj(x_prev), "\n")
   if (abs(obj(x_cur) - obj(x_prev)) < eps) {
    #cat("Converges at", k, "iteration.  \n")
    #cat("objective function is ", obj(x_cur), "\n")
    break; 
  }
  if (k == max.iter) {
    #cat("max number of iteration reached! \n")
    #cat("objective function is ", obj(x_cur), "\n")
  }
}
x_cur
}



prox_grad_nonconvex <- function(obj_g, obj, grad_g, prox_h, x0=NULL, dim=NULL, max.iter = 2e3, beta = .8, t = 1e-2,eps = 1e-5){
  if (is.null(x0)) {
    if (is.null(dim)) stop("Please specify either an initial solution or input dimension")
      x0 = rep(0,dim)
  }
  else {
    if (is.null(dim)) { dim = length(x0) }
    else { if (dim != length(x0)) stop("dimension of initial solution must have dimension dim") }
  }
  x_prev = x_cur = z = x0
  delta_x = rep(0,dim)
  for (k in 1:max.iter) {
    x_prev = x_cur
    x_cur = prox_h(x_prev-t*grad_g(x_prev), t)
    delta_x = x_cur - x_prev
    # cat("step size is", t, "\n")
    # cat("objective function is ", obj(x_prev), "\n")
    while (obj_g(x_cur) > (obj_g(x_prev) + crossprod(grad_g(x_prev), delta_x) + sum(delta_x^2) / (2*t))) {
      t = t * beta
      x_cur = prox_h(x_prev-t*grad_g(x_prev), t)
      delta_x = x_cur - x_prev
    }
    ####### check convergence #######
    if (abs(obj(x_cur) - obj(x_prev)) < eps) {
      #cat("Converges at", k, "iteration.  \n")
      #cat("objective function is ", obj(x_cur), "\n")
      break; 
    }
    if (k == max.iter) {
      #cat("max number of iteration reached! \n")
      #cat("objective function is ", obj(x_cur), "\n")
    }
  }
  
  x_cur
}







