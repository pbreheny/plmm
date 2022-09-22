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
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' my_fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' LL <- logLik_nonnull(fit = my_fit, SUX = my_fit$SUX, SUy = my_fit$SUy, 
#' S = my_fit$S, eta = my_fit$eta)
#' head(LL)
#' # See pp. 16-18 in A.R.'s thesis for details 

logLik_nonnull <- function(fit, SUX, SUy, S, eta){

  n <- dim(SUX)[1]
  beta_vals <- coef.plmm(fit, drop = FALSE) # don't drop - keep names
  n_betas <- dim(coef.plmm(fit))[1]
  # choose only beta values from nonsingular SNPs (predictors)
  # NB: the intercept is included here 
  b <- beta_vals[fit$ns_idx,] 

  # evaluate log determinant
  Sd <- eta * S + (1 - eta)
  ldet <- sum(log(Sd))

  # evaluate log likelihood at each value of lambda
  ll <- apply(b,
              2,
              function(bb){
    ss <- (1/n) * crossprod(SUy - SUX %*% bb)
    -0.5 * (n*log(2*pi) + ldet + n + n*log(ss))
    }
  )
  
  # pad NAs where there are monomorphic values 
  ll_padded <- rep(NA_integer_, n_betas)
  names(ll_padded) <- names(beta_vals)
  ll_padded[fit$ns_idx] <- ll

  return(ll_padded)

}

