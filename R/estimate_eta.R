#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @export
#' @keywords internal
estimate_eta <- function(s, U, y){

  # optimize over \eta and \beta_0 
  opt <- optim(par = c(0.2, 0.1),
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
#' @param s the vector of singular values of data X 
#' 
#' @keywords internal
null_model <- function(params, y, U, s){
  
  # name parameters
  eta <- params[1]
  beta0 <- params[2]
  
  # create intercept (for null model, this is the only predictor)
  n <- length(y)
  intcpt <- rep(1, n)
  
  # get \Sigma^2_{-1/2} piece
  w2 <- (eta*s) + (1 - eta)
  w <- w2^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  rot_y <- wUt %*% y 
  rot_intcpt <- wUt %*% intcpt
  rot_mean <- drop(rot_intcpt*beta0)
  
  # distribution of null model on rotated scale
  res <- mvtnorm::dmvnorm(x = drop(rot_y),
               mean = rot_mean,
               log = TRUE) # use log scale for numerical stability
  
  # -0.5*(n*log(2*pi) + sum(log())) # thinking about writing out the density myself...

  return(-1*res) # optim() does minimization, so I need the neg. log. lik.

}


#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param beta0 The coefficient corresponding to the intercept (in null model, 
#' this is the only coefficient).
#' @param K Matrix to use as relatedness matrix.
#' @keywords internal
test_eta_estimation <- function(sig_s, sig_eps, beta0, K){
  
  # true_eta <- sig_s/(sig_s + sig_eps)
  intcpt <- rep(1, nrow(K))
  
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  eps <- mvtnorm::rmvnorm(n = 1,
                          sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  
  y <- intcpt*beta0 + u + eps # null model = intercept only model s

  eig_K <- eigen(K)
  tmp <- estimate_eta(s = eig_K$values, U = eig_K$vectors, y = y)
  
  return(list(hat_eta = tmp$eta,
              hat_beta0 = tmp$beta0))
}


# pick up here: alternative setup -- write a nLL function and call this from estimate_eta

