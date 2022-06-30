#' Evaluate the negative log-likelihood of a null Gaussian penalizedLMM model
#'
#' This function allows you to evaluate the negtive log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param eta The proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR. Sometimes referred to as the narrow-sense heritability.
#' @param Uy The the continuous outcome, y, rotated by the eigenvectors of the similarity matrix, K.
#' @param S The eigenvalues of the similarity matrix, K.
#' @export
#' 
#' @examples 
#' admix$K <- admix$X%*%t(admix$X) # create an estimated covariance matrix 
#' # NB: this is an estimate of K 
#' ev <- eigen(admix$V)
#' U <- ev$vectors
#' fit <- plmm(X = admix$X, y = admix$y, V = admix$V, penalty = "MCP")
#' (logLik(eta = fit$eta, Uy = U%*%admix$y, S = ev$values ))

logLik <- function(eta, Uy, S){

  n <- dim(Uy)[1]

  # evaluate log determinant
  Sd <- eta * S + (1 - eta)
  ldet <- sum(log(Sd))  # log of product = sum of the logs

  # evaluate the variance
  Sdi <- 1/Sd
  Uy <- as.vector(Uy)
  ss <- (1/n) * sum(Uy*Uy*Sdi)

  # evalue the negative log likelihood
  nLL <- 0.5*(n*log(2*pi) + ldet + n + n*log(ss))

  return(nLL)

}
