#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @param init.vals Vector of length 2 with the intial values for eta and beta0, in that order. Defaults to c(0.2, 0.1). Passed as `par` to `optim()`.
#' @export
#' @keywords internal
estimate_eta <- function(s, U, y, init.vals = c(0.2, 0.1)){

  # optimize over \eta and \beta_0 
  opt <- optim(par = init.vals,
                  fn = null_model,
                  y = y,
                  s = s,
                  U = U,
                  method = "L-BFGS-B",
               # use box constraints to ensure 0 < \hat \eta < 1
               lower = c(1e-5, -Inf),
               upper = c(0.99999, Inf))

  return(list(eta = opt$par[1],
              beta0 = opt$par[2]))
}

#' a helper function for 2-dim optimization
#' @param params A vector of 3 elements, corresponding to eta and beta_0 
#' The latter is the coefficient of the null model
#' @param y A vector of outcomes 
#' @param U The left singular vectors of data X 
#' @param s the vector of singular values of data K 
#' 
#' @keywords internal
null_model <- function(params, y, U, s){
  
  # name parameters
  eta <- params[1]
  beta0 <- params[2]

  # create intercept (for null model, this is the only predictor)
  n <- length(y)
  intcpt <- rep(1, n)
  
  # get \Sigma^2_{-1/2} piece - call it wUt (for \W \U^\top)
  w2 <- ((eta*s) + (1 - eta))
  w <- w2^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  rot_y <- wUt %*% y
  rot_intcpt <- (wUt %*% intcpt)
  rot_mean <- drop(rot_intcpt*beta0)
  
  # distribution of null model on rotated scale
  res <- mvtnorm::dmvnorm(x = drop(rot_y),
               mean = rot_mean,
               log = TRUE) # use log scale for numerical stability
  return(-1*res) # optim() does minimization, so I need the neg. log. lik.

  # working out the nLL myself: 
  # constant <- n*log(2*pi)
  # rot_resid <- crossprod(rot_y - rot_mean)
  # nll <- 0.5*(constant + rot_resid)
  # return(nll)

}


#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param beta0 The coefficient corresponding to the intercept (in null model, 
#' this is the only coefficient).
#' @param K Matrix to use as relatedness matrix.
#' @param ... Additional args to pass into `estimate_eta()`
#' @keywords internal
test_eta_estimation <- function(sig_s, sig_eps, beta0, K, ...){
  
  # true_eta <- sig_s/(sig_s + sig_eps)
  intcpt <- rep(1, nrow(K))
  
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  eps <- mvtnorm::rmvnorm(n = 1,
                          sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  
  y <- intcpt*beta0 + u + eps # null model = intercept only model s

  eig_K <- eigen(K)
  tmp <- estimate_eta(s = eig_K$values, U = eig_K$vectors, y = y, ...)
  
  return(list(hat_eta = tmp$eta,
              hat_beta0 = tmp$beta0))
}