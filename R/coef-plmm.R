#' Coef method for "plmm" class
#'
#' @param object An object of class "plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param drop Logical.
#' @param ... Additional arguments.
#' @export
#' 
#' @examples 
#' fit <- plmm(admix$X, admix$y)
#' (coef.plmm(fit)[1:10, 1:5])

coef.plmm <- function(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...){
  if (!missing(lambda)) {
    if (max(lambda) > max(object$lambda) | min(lambda) <
        min(object$lambda)) {
      stop("Supplied lambda value(s) are outside the range of the model fit.",
           call. = FALSE)
    }
    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind%%1
    beta_vals <- (1 - w) * object$beta_vals[, l, drop = FALSE] + w * object$beta_vals[, r, drop = FALSE]
    colnames(beta_vals) <- lamNames(lambda)
  }
  else beta_vals <- object$beta_vals[, which, drop = FALSE]
  if (drop){
    return(drop(beta_vals))
  } else{
    return(beta_vals)
  }
}
