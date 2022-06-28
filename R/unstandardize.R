#' Unstandardize coefficient values
#'
#' This function allows you to unstandardize coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Original, non-standardized design matrix without an intercept column. Necessary for properly defining the dimensions of beta in cases where singular columns are present.
#' @param std_X Standardized design matrix. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export
#' 
#' @examples 
#' admix$K <- (admix$X%*%t(admix$X))/ncol(admix$X) # create an estimated covariance matrix 
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' # TODO: this throws an error: new_betas <- unstandardize(b = fit$beta, X = admix$X, std_X = scale(admix$X))

unstandardize <- function(b, X, std_X, intercept = TRUE) {
  ns <- attr(std_X, 'nonsingular')
  center <- attr(std_X, 'center')[ns]
  scale <- attr(std_X, 'scale')[ns]
  if (intercept){
    beta <- matrix(0, nrow = ncol(X) + 1, ncol = ncol(b))
    a <- b[1, , drop = FALSE]
    b <- b[-1, , drop=FALSE]
    bb <- b / scale
    beta[1 + ns,] <- bb
    beta[1,] <- a - crossprod(center, bb)
  } else {
    beta <- matrix(0, nrow = ncol(X), ncol = ncol(b))
    bb <- b / scale
    beta[ns, ] <- bb
  }
  return(beta)
}
