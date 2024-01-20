#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @export
#' @keywords internal
estimate_eta_null_mean_0 <- function(s, U, y){
  
  # optimize over \eta 
  opt <- optimize(f = null_model_mean_0,
                  interval = c(1e-5, 0.99999),
                  y = y,
                  s = s,
                  U = U)
  
  return(list(eta = opt$minimum,
              objective = opt$objective))
}

#' a helper function for 2-dim optimization
#' @param params A vector of 3 elements, corresponding to eta and beta_0 
#' The latter is the coefficient of the null model
#' @param y A vector of outcomes 
#' @param U The left singular vectors of data X 
#' @param s the vector of singular values of data K 
#' 
#' @keywords internal
null_model_mean_0 <- function(eta, y, U, s){
  # get \Sigma^2_{-1/2} piece - call it wUt (for \W \U^\top)
  w2 <- ((eta*s) + (1 - eta))
  w <- w2^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  rot_y <- wUt %*% y
  
  # distribution of null model on rotated scale
  res <- mvtnorm::dmvnorm(x = drop(rot_y),
                          # use default mean and sd; we have a standard normal here
                          log = TRUE # use log scale for numerical stability
  )
  return(-1*res) # optim() does minimization, so I need the neg. log. lik.
  
}


#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param beta0 The coefficient corresponding to the intercept (in null model, 
#' this is the only coefficient).
#' @param K Matrix to use as relatedness matrix.
#' @param ... Additional args to pass into `estimate_eta()`
#' @keywords internal
test_eta_estimation_null_mean_0 <- function(sig_s, sig_eps, beta0, K, ...){
  
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  eps <- mvtnorm::rmvnorm(n = 1,
                          sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  
  y <- u + eps # null model here = no predictors, not even intercept
  
  eig_K <- eigen(K)
  tmp <- estimate_eta_null_mean_0(s = eig_K$values, U = eig_K$vectors, y = y, ...)
  
  return(list(hat_eta = tmp$eta,
              objective = tmp$objective))
}