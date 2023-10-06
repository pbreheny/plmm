#' Evaluate the negative log-likelihood of a null Gaussian penalizedLMM model
#'
#' This function allows you to evaluate the negtive log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param eta The proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR. Sometimes referred to as the narrow-sense heritability.
#' @param rot_y The the continuous outcome, y, rotated by the eigenvectors of the similarity matrix, K.
#' @param s The eigenvalues of the similarity matrix, K.
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' ev <- eigen(admix$K)
#' U <- ev$vectors
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' (log_lik(eta = fit$eta, rot_y = U%*%admix$y, s = ev$values ))
#' }

log_lik <- function(eta, rot_y, s){

  n <- dim(rot_y)[1]

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
