#' Predict method for "plmm" class
#'
#' @param object An object of class "plmm."
#' @param X Design matrix used for computing predicted values if requested.
#' @param type a character argument indicating what should be returned.
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param ... Additional arguments.
#' @export

## from ncvreg
predict.plmm <- function(object, X, type=c("response", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), ...) {
  ### I think this might be wrong/need to be more involved than simply X %*% beta-hat bc of ranef
  ### Need to use BLUPs blerrrrg

  # need better way to deal with intercept...
  type <- match.arg(type)
  beta <- coef.plmm(object, lambda=lambda, which=which, drop=FALSE) # includes intercept
  if (type=="coefficients") return(beta)
  # alpha <- beta[1, , drop = FALSE]
  # beta <- beta[-1, , drop=FALSE]
  if (type=="nvars") return(apply(beta[-1, , drop=FALSE] !=0, 2, sum)) # don't count intercept
  if (type=="vars") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  eta <- X %*% beta
  if (type=="response") return(drop(eta))
}
