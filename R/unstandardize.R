#' Unstandardize coefficient values
#'
#' This function allows you to unstandardize coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Design matrix which *includes* the intercept column if present. May include clinical covariates and other non-SNP data. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export

unstandardize <- function(b, X, intercept = TRUE) {
  beta <- matrix(0, nrow = nrow(b), ncol = ncol(b))
  if (intercept){
    ns <- attr(X, 'nonsingular')
    a <- b[1, , drop = FALSE]
    b <- b[-1, , drop=FALSE]
    bb <- b / attr(X, 'scale')[ns]
    beta[1 + ns,] <- bb
    beta[1,] <- a - crossprod(attr(X, 'center')[ns], bb)
  } else {
    ns <- attr(X, 'nonsingular')
    center <- attr(X, 'center')[ns]
    scale <- attr(X, 'scale')[ns]
    bb <- b / attr(X, 'scale')[ns]
    beta[ns, ] <- bb
  }
  return(beta)
}
