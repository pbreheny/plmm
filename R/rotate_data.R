#' Compute transformed X and y for subsequent plmm rescaling and fitting
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param K Estimated or known similarity matrix. By default, K is the realized relationship matrix, \eqn{\frac{1}{p}XX^T}, where \eqn{p} is the number of columns in X
#' @param intercept Logical flag for whether an intercept should be included.
#' @param rotation Logical flag to indicate whether the rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @param eta_star Optional arg to input a specific eta term rather than estimate it from the data. If v is a known matrix, this should be 1.
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' std_X <- scale(admix$X)
#' K <- relatedness_mat(std_X)
#' rotated_dat <- rotate_data(std_X, admix$y, K)

rotate_data <- function(X, y, K = NULL, intercept = TRUE, rotation = TRUE, eta_star){
  # Coersion
  S <- U <- UU <- eta <- NULL

  # Calculate RRM
  if(is.null(K)){
    p <- ncol(X)
    K <- relatedness_mat(X)
  }
  
  ## Calculate U, S, eta
  if (rotation){
    if (missing(eta_star)){
      c(S, U, eta) %<-% plmm_null(K = K, y = y) # estimate eta if needed 
    } else {
      c(S, U) %<-% svd(K)[1:2]
      eta <- eta_star
    }
    W <- diag((eta * S + (1 - eta))^(-1/2))
  } else {
    # no rotation option for testing
    S <- rep(1, nrow(X))
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
