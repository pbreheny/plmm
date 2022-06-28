#' Evaluate the negative log-likelihood of a non-null Gaussian plmm model
#'
#' @param fit An object of class plmm_fit.
#' @param SUX Rotated design matrix including rotated intercept and unpenalized columns, if present.
#' @param SUy Rotated outcome vector.
#' @param S Eigenvalues from similarity matrix used for model fitting.
#' @param eta Estimated $eta$ value from object fit.
#' @export
#' 
#' @examples 
#' admix$K <- (1/nrow(admix$X))*tcrossprod(admix$X, admix$X) # create an estimated covariance matrix 
#' # NB: this is an estimate of K 
#' K_svd <- svd(K)
#' U <- K_svd$u
#' D <- K_svd$d
#' # WORK IN PROGRESS #FIXME need to iron out this example
#' # See pp. 16-18 in A.R.'s thesis 
#' ev <- eigen(admix$V)
#' U <- ev$vectors
#' S <- ev$values
#' UX <- U%*%admix$X
#' Uy <- U%*%admix$y
#' fit <- plmm(X = admix$X, y = admix$y, V = admix$V, penalty = "MCP")
#' (plmm_nll_nonnull(fit = my_fit, SUX = S%*%UX, SUy = S%*%admix$y, S = S, eta = my_fit$eta))

plmm_nll_nonnull <- function(fit, SUX, SUy, S, eta){

  n <- dim(SUX)[1]
  beta <- coef.plmm(fit)

  # evaluate log determinant
  Sd <- eta * S + (1 - eta)
  ldet <- sum(log(Sd))

  ll <- apply(beta, 2, function(bb){
    ss <- (1/n) * crossprod(SUy - SUX %*% bb)
    -0.5 * (n*log(2*pi) + ldet + n + n*log(ss))
  })

  return(ll)

}

