#' Coef method for a list used in cross-validation (e.g. the \code{fit.i} within \code{cvf})
#'
#' @param fit A list with the components returned by `plmm_fit`
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param drop Logical.
#' @param ... Additional arguments.
#' @keywords internal
#' 

coef.list <- function(fit, lambda, which, drop, ...){
  # error check for supplied lambda value 
  if (!missing(lambda)) {
    if (max(lambda) > max(fit$lambda) | min(lambda) <
        min(fit$lambda)) {
      stop("Supplied lambda value(s) are outside the range of the model fit.",
           call. = FALSE)
    }
    ind <- stats::approx(fit$lambda, seq(fit$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind%%1
    beta_vals <- (1 - w) * fit$b[, l, drop = FALSE] + w * fit$b[, r, drop = FALSE]
    
    # format dim. names
    if(is.null(dim(beta_vals))) {
      # case 1: beta_vals is a vector 
      names(beta_vals) <- lamNames(lambda)
    } else {
      # case 2: beta_vals is a matrix
      colnames(beta_vals) <- lamNames(lambda)
    }
    
  }
  else beta_vals <- fit$b[, which, drop = FALSE]
  
  if (drop){
    return(drop(beta_vals))
  } else{
    return(beta_vals)
  }
}