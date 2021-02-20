#' Predict method for plmm class
#'
#' @param object An object of class plmm.
#' @param newX Design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what should be returned.
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param XX Optional argument. Original design matrix (not including intercept column) from object. Required if \code{type == 'individual'}.
#' @param y Optional argument. Original continuous outcome vector from object. Required if \code{type == 'individual'}.
#' @param U Optional argument. Eigenvectors from the similarity matrix from object. Required if \code{type == 'individual'}.
#' @param S Optional argument. Eigenvalues from the similarity matrix from object. Required if \code{type == 'individual'}.
#' @param eta Optional argument. Estimated $eta$ value from object. Required if \code{type == 'individual'}.
#' @param covariance Optional argument. $q times n$ covariance matrix between new and old observations. Required if \code{type == 'individual'}.
#' @param ... Additional arguments.
#' @export

predict.plmm <- function(object, newX, type=c("response", "individual", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), XX, y, U, S, eta, covariance, ...) {
  type <- match.arg(type)
  beta <- coef.plmm(object, lambda=lambda, which=which, drop=FALSE) # includes intercept
  if (type=="coefficients") return(beta)
  if (type=="nvars") return(apply(beta[-1, , drop=FALSE] !=0, 2, sum)) # don't count intercept
  if (type=="vars") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  Xbeta <- cbind(1, newX) %*% beta
  if (type=="response") return(drop(Xbeta))
  if (type == "individual"){
   if (missing(XX) || missing(y) || missing(U) || missing(S) || missing(eta) || missing(covariance)){
     stop('Original XX, y, U, S, eta, and covariance must be supplied if type is individual')
   }
   # can't just use the rotated y and x here - need to scale by inverse of V, not sqrt(V)
   ranef <- covariance %*% U %*% diag((1 + eta * (S - 1))^(-1)) %*% t(U) %*% (y - cbind(1, XX) %*% beta)
   # print(eta)
   blup <- Xbeta + ranef
   return(blup)
  }
}
