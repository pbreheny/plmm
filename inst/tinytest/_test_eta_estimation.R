#' a function to write tests of estimate_eta
#' @param sig_s Variance attributable to structure 
#' @param sig_eps Variance of random error 
#' @param K Matrix to use as relatedness matrix.
#' @param ... Additional args to pass into `estimate_eta()`
#' @keywords internal
test_eta_estimation <- function(sig_s, sig_eps, K, ...){
  
  # Note: true_eta <- sig_s/(sig_s + sig_eps)
  
  # simulate data
  intcpt <- rep(1, nrow(K))
  
  u <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_s*K) |> drop()
  
  eps <- mvtnorm::rmvnorm(n = 1,
                          sigma = sig_eps*diag(nrow = nrow(K))) |> drop()
  
  y <- intcpt + u + eps # null model = intercept only model s
  
  eig_K <- eigen(K)
  
  # check signs 
  # nz <- which(eig_K$values > 0.00000001)
  # sign_check <- flip_signs(X = K,
  #                          U = eig_K$vectors[,nz],
  #                          V = eig_K$vectors[,nz], 
  #                          d = eig_K$values[nz])
  # U <- sign_check$U
  
  # estimate eta
  tmp <- estimate_eta(s = eig_K$values, U = eig_K$vectors, y = y, ...)
  
  return(tmp)
}


K <- relatedness_mat(admix$X)
hat_eta <- rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 2,
                             sig_eps = 1,
                             K = K)
  hat_eta[i] <- res
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta) 
