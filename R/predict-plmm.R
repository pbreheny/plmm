#' Predict method for "plmm" class
#'
#' @param object An object of class "plmm."
#' @param newX Design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what should be returned.
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param X Optional argument. Original design matrix (not including intercept column) from object fit. Required if \code{type == 'individual'}.
#' @param y Optional argument. Original continuous outcome vector from object fit. Required if \code{type == 'individual'}.
#' @param eigenvectors Optional argument. Eigenvectors from the similarity matrix from object fit.  Required if \code{type == 'individual'}.
#' @param eigenvalues Optional argument. Eigenvalues from the similarity matrix from object fit. Required if \code{type == 'individual'}.
#' @param eta Optional argument. Estimated $eta$ value from object fit. Required if \code{type == 'individual'}.
#' @param covariance Optional argument. $q times n$ covariance matrix between new and old observations. Required if \code{type == 'individual'}.
#' @param ... Additional arguments.
#' @export

## from ncvreg
predict.plmm <- function(object, newX, type=c("response", "individual", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), X, y, eigenvectors, eigenvalues, eta, covariance, ...) {
  type <- match.arg(type)
  beta <- coef.plmm(object, lambda=lambda, which=which, drop=FALSE) # includes intercept
  if (type=="coefficients") return(beta)
  if (type=="nvars") return(apply(beta[-1, , drop=FALSE] !=0, 2, sum)) # don't count intercept
  if (type=="vars") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  Xbeta <- cbind(1, newX) %*% beta
  if (type=="response") return(drop(Xbeta))
  if (type == "individual"){
   if (missing(X) || missing(y) || missing(eigenvectors) || missing(eigenvalues) || missing(eta) || missing(covariance)){
     stop('Original X, y, eigenvectors, eigenvalues, eta, and covariance must be supplied if type is individual')
   }
   # can't just use the rotated y and x here - need to scale by inverse of V, not sqrt(V)
   # ranef <- covariance %*% solve(eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors) * eta + diag(length(y)) * (1 - eta)) %*% (y - X %*% beta)
   # ranef <- covariance %*% eigenvectors %*% solve(diag(eigenvalues) * eta + diag(length(y)) * (1 - eta)) %*% t(eigenvectors) %*% (y - X %*% beta)
   # ranef <- covariance %*% eigenvectors %*% diag((eigenvalues * eta + 1 - eta)^(-1)) %*% t(eigenvectors) %*% (y - X %*% beta)
   # make sure this vector is the correct length if singular
   eigenvalues <- c(eigenvalues, rep(0, nrow(X) - length(eigenvalues)))
   ranef <- covariance %*% eigenvectors %*% diag((1 + eta * (eigenvalues - 1))^(-1)) %*% t(eigenvectors) %*% (y - cbind(1, X) %*% beta)

   # print(eta)
   blup <- Xbeta + ranef
   return(blup)
  }
}
