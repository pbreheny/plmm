#' Unscale coefficient values
#'
#' This function allows you to unscale coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Unscaled, rotated design matrix *with* an intercept column if present. Necessary for properly defining the dimensions of beta in cases where singular columns are present.
#' @param XX Scaled, rotated design matrix. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export

unscale <- function(b, X, XX, intercept = TRUE) {
  beta <- matrix(0, nrow = ncol(X), ncol = ncol(b))
  ns <- attr(XX, 'nonsingular')
  scale <- attr(XX, 'scale')[ns]
  if (intercept){
    beta[1,] <- b[1, , drop = FALSE]
    b <- b[-1, , drop=FALSE]
    bb <- b / scale
    beta[1 + ns,] <- bb
  } else {
    bb <- b / scale
    beta[ns, ] <- bb
  }
  return(beta)
}
