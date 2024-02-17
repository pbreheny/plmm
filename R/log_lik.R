#' Evaluate the negative log-likelihood of an intercept-only Gaussian plmm model
#'
#' This function allows you to evaluate the negative log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param eta The proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR. Sometimes referred to as the narrow-sense heritability.
#' @param n The number of observations 
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @param rot_y Optional: if y has already been rotated, then this can be supplied. This option is designed for the `log_lik()` call within `gic()`, where a `fit` object is being supplied.
#' @keywords internal
#' 

log_lik <- function(eta, n, s, U, y, rot_y = NULL){

  # first, the constant (comes from 1st term in derivation)
  constant <- n*log(2*pi)
  
  # TODO: think about simplifying what follows with `rotate()`
  # we will need the sum of the nonzero values from the diagonal matrix of weights
  w2 <- ((eta*s) + (1 - eta))
  
  # get w2 on the log scale 
  sum_det_log <- sum(log(w2))
  
  # rotate y 
  w <- w2^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  
  
  if(is.null(rot_y) & !(missing(y))){
    rot_y <- wUt %*% y
  } 
  
  
  # rotate the intercept (this is the only term in the null model)
  intcpt <- rep(1, n)
  rot_intcpt <- (wUt %*% intcpt)
  
  # calculate the term representing the \hat\beta(\eta) MLE
  intcpt_crossprod <- crossprod(rot_intcpt)
  intcpt_y_crossprod <- crossprod(rot_intcpt, rot_y)
  hat_beta_mle <- drop(intcpt_y_crossprod/intcpt_crossprod)
  
  # using hat_beta_mle, calculate the quadratic term from the log likelihood
  rot_e <- rot_y - (rot_intcpt*hat_beta_mle) # e = 'error', y - mean
  rot_sqe <- crossprod(rot_e) # mse = sq. error
  quad_term <- (1/n)*rot_sqe/sum(w2)
  
  # put all the pieces together -- evaluate the **negative** log likelihood
  # NB: keep constant here to be consistent with log_lik.lm() method
  nLL <- 0.5*(constant + n*log(quad_term) + sum_det_log + n)
  
  return(drop(nLL))

}
