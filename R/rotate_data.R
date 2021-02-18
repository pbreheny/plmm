#' Compute transformed X and y for subsequent plmm rescaling and fitting
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param intercept Logical flag for whether an intercept should be included.
#' @param rotation Logical flag to indicate whether the rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @param eta_centerY Optional logical flag to indicate whether eta should be estimated using a centerd outcome variable, y. Defaults to FALSE.
#' @param eta_star Optional arg to input a specific eta term rather than estimate it from the data.
#' @param V Optional arg to input a specified covariance structure for y, used to rotate the data, into the model.
#' @importFrom zeallot %<-%
#' @export

rotate_data <- function(X, y, X_for_K, intercept, rotation = TRUE, eta_centerY = FALSE, eta_star = NULL, V = NULL){
  # Coersion
  d <- U <- UU <- eta <- NULL

  ## Calculate U, d, eta
  if (rotation){
    if (!is.null(V) && dim(V)[1] == nrow(X)){
      c(d, U, UU) %<-% svd(V)
      W <- diag(1/sqrt(d))
      eta <- NA
    } else if (is.null(V)){
      if (eta_centerY){
        c(d, U, eta) %<-% plmm_null(X_for_K, y - mean(y))
      } else{
        c(d, U, eta) %<-% plmm_null(X_for_K, y)
      }
      # still compute U and d but override eta-hat with eta_star if supplied
      if (!is.null(eta_star)) eta <- eta_star
      W <- diag((eta * d + (1 - eta))^(-1/2))
    }
  } else {
    # no rotation option for testing
    d <- diag(nrow(X))
    U <- diag(nrow(X))
    W <- diag(nrow(X))
    eta <- NA
  }

  ## Rotate data
  if (intercept){
    SUX <- W %*% crossprod(U, cbind(1, X))
  } else {
    SUX <- W %*% crossprod(U, X)
  }
  SUy <- drop(W %*% crossprod(U, y))
  val <- list(X = SUX,
              y = SUy,
              eta = eta,
              U = U,
              d = d)
  return(val)
}
