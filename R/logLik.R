#' Evaluate the negative log-likelihood of a null Gaussian penalizedLMM model
#'
#' This function allows you to evaluate the negtive log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param eta The proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR. Sometimes referred to as the narrow-sense heritability.
#' @param Uy The the continuous outcome, y, rotated by the eigenvectors of the similarity matrix, K.
#' @param S The eigenvalues of the similarity matrix, K.
#' @export
#' 
#' @examples 
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' ev <- eigen(admix$K)
#' U <- ev$vectors
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' (logLik(eta = fit$eta, Uy = U%*%admix$y, S = ev$values ))

logLik <- function(eta, Uy, S){
  
  n <- dim(Uy)[1]
  
  # evaluate log determinant
  Sd <- eta * S + (1 - eta) # Thesis page 18 W^(-2)
  ldet <- sum(log(Sd))  # log of product = sum of the logs
  
  # evaluate the variance
  Sdi <- 1/Sd
  Uy <- as.vector(Uy)
  ss <- sum(Uy*Uy*Sdi)
  
  # evalue the negative log likelihood
  nLL <- 0.5*(ldet + log(ss)) # ss is on log scale for numerical stability reason 
  
  return(nLL)
  
}
