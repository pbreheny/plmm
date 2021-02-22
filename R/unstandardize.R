#' Unstandardize coefficient values
#'
#' This function allows you to unstandardize coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Original, non-standardized design matrix without an intercept column. Necessary for properly defining the dimensions of beta in cases where singular columns are present.
#' @param XX Standardized design matrix. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export

unstandardize <- function(b, X, XX, intercept = TRUE) {
  ns <- attr(XX, 'nonsingular')
  center <- attr(XX, 'center')[ns]
  scale <- attr(XX, 'scale')[ns]
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
