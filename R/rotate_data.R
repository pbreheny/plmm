#' Compute transformed X and y for subsequent plmm rescaling and fitting
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param V Similarity matrix. This may be a known value supplied to plmm as V, or estimated based on X_for_K
#' @param eta_star Optional arg to input a specific eta term rather than estimate it from the data. If v is a known matrix, this should be 1.
#' @param intercept Logical flag for whether an intercept should be included.
#' @param rotation Logical flag to indicate whether the rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @importFrom zeallot %<-%
#' @export

rotate_data <- function(X, y, V, eta_star, intercept, rotation = TRUE){
  # Coersion
  S <- U <- UU <- eta <- NULL

  ## Calculate U, S, eta
  if (rotation){
    if (missing(eta_star)){
      c(S, U, eta) %<-% plmm_null(V, y)
    } else {
      c(S, U, UU) %<-% svd(V)
      eta <- eta_star
    }
    W <- diag((eta * S + (1 - eta))^(-1/2))
  } else {
    # no rotation option for testing
    S <- diag(nrow(X))
    U <- diag(nrow(X))
    W <- diag(nrow(X))
    eta <- 1
  }

  ## Rotate data
  if (intercept){
    SUX <- W %*% crossprod(U, cbind(1, X))
  } else {
    SUX <- W %*% crossprod(U, X)
  }
  SUy <- drop(W %*% crossprod(U, y))
  val <- list(SUX = SUX,
              SUy = SUy,
              eta = eta,
              U = U,
              S = S)
  return(val)
}
