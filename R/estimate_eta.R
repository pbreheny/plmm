#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param n The number of observations 
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' 
#' @keywords internal
estimate_eta <- function(n, s, U, y, eta_star){

  opt <- stats::optimize(f=log_lik,
                         c(0.01, 0.99),
                         n = n,
                         s = s,
                         U = U,
                         y = y)
  
  eta <- opt$minimum 
  
  
  return(eta)
}

#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param K Matrix to use as relatedness matrix.
#' @param intercept Logical: should the outcome be simulated with an intercept? Defaults to TRUE. 
#' @param center_y  Logical: should y be centered? Defaults to FALSE. (Note: this option only makes sense if intercept = TRUE)
#' @param y_skew Numeric value: if nonzero, this will control how much skewness will be added to the noise. Defaults to 0. 
#' @param return_y Logical: should simulated outcome values be returned? Defaults to FALSE.
#' @param ... Additional args to pass into `estimate_eta()`
#' @keywords internal
test_eta_estimation <- function(sig_s, sig_eps, K, intercept = TRUE,
                                center_y = FALSE, y_skew = 0, return_y = FALSE, ...){
  
  # Note: true_eta <- sig_s/(sig_s + sig_eps)
  
  # simulate data
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  
  if (!(y_skew == 0)) {
    t <- runif(1)
    if (t <= 0.1) {
      eps <- mvtnorm::rmvnorm(n = 1,
                              sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
      eps <- eps + y_skew
    }
  } else {
    eps <- mvtnorm::rmvnorm(n = 1,
                            sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  }
  
  y <- u + eps
  if (intercept){
    intcpt <- rep(1, nrow(K))
    y <- y + intcpt
    
    if (center_y){
      y <- y - mean(y)
    }
  }

  eig_K <- eigen(K)
  
  # check signs 
  # TODO: determine if we need this here
  # nz <- which(eig_K$values > 0.00000001
  # sign_check <- flip_signs(X = K,
  #                          U = eig_K$vectors[,nz],
  #                          V = eig_K$vectors[,nz], 
  #                          d = eig_K$values[nz])
  # U <- sign_check$U
  
  # estimate eta
  eta <- estimate_eta(n = length(y),
                      s = eig_K$values,
                      U = eig_K$vectors,
                      y = y,
                       ...)
  if (return_y){
    return(structure(list(y = y, eta = eta)))
  } else {
    return(eta)
  }
  
}