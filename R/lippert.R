#' the estimation of eta using the Lippert (2011) derivation
lippert_estimate_eta <- function(n, s, U, y, eta_star){
  
  # coercion
  eta <- NULL
  
  # estimate eta 
  rot_y <- crossprod(U, y)
  opt <- stats::optimize(f=log_lik, c(0.01, 0.99), rot_y=rot_y, s=s, n=nrow(U))
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
  # TODO: double check the derivation for this. Do we need the factor of n 
  #   in the last term? 
  return(nLL)
  
}
