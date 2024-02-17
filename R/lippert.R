#' the estimation of eta using the Lippert (2011) derivation
#'
#' @param n Number of observations 
#' @param s Eigenvalues of K
#' @param U Eigenvectors of K
#' @param y Outcome vector 
lippert_estimate_eta <- function(n, s, U, y){
  
  # coercion
  eta <- NULL
  
  # estimate eta 
  rot_y <- crossprod(U, y)
  opt <- stats::optimize(f=lippert_loglik, c(0.01, 0.99), rot_y=rot_y, s=s, n=nrow(U))
  eta <- opt$minimum 
  
  return(eta)
}

#' The optimization function from Lippert (2011)
#' 
lippert_loglik <- function(eta, rot_y, s, n){
  
  # evaluate log determinant
  sd <- eta * s + (1 - eta)
  ldet <- sum(log(sd))  # log of product = sum of the logs
  
  # evaluate the variance
  sdi <- 1/sd 
  rot_y <- as.vector(rot_y)
  ss <- (1/n) * sum(rot_y*rot_y*sdi)
  
  # evaluate the negative log likelihood
  # NB: keep constant here to be consistent with log_lik.lm() method
  nLL <- 0.5*(n*log(2*pi) + ldet + n + n*log(ss))

  return(nLL)
  
}

#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param K Matrix to use as relatedness matrix.
#' @param ... Additional args to pass into `estimate_eta()`
#' @keywords internal
lippert_test_eta_estimation <- function(sig_s, sig_eps, K, ...){
  
  # Note: true_eta <- sig_s/(sig_s + sig_eps)
  
  # simulate data
  intcpt <- rep(1, nrow(K))
  
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  eps <- mvtnorm::rmvnorm(n = 1,
                          sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  
  y <- intcpt + u + eps # null model = intercept only model s
  
  eig_K <- eigen(K)
  
  # estimate eta
  eta <- lippert_estimate_eta(n = length(y),
                      s = eig_K$values,
                      U = eig_K$vectors,
                      y = y)
  
  return(eta)
}
