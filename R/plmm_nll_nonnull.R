#' Evaluate the negative log-likelihood of a non-null Gaussian plmm model
#'
#' @param fit An object of class plmm_fit.
#' @param SUX Rotated design matrix including rotated intercept and unpenalized columns, if present.
#' @param SUy Rotated outcome vector.
#' @param S Eigenvalues from similarity matrix used for model fitting.
#' @param eta Estimated $eta$ value from object fit.
#' @export

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

