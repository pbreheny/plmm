#' Coef method for "plmm" class
#'
#' @param object An object of class "plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param drop Logical.
#' @param ... Additional arguments.
#' 
#' @rdname coef.plmm
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' fit <- plmm(admix$X, admix$y)
#' (coef.plmm(fit)[1:10, 1:5])
#' # TODO: in R CMD CHECK, this throws the error: 
#' # could not find function "coef.plmm"
#' }
#' 
#' 

coef.plmm <- function(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...){
  # error check for supplied lambda value 
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
    
    # format dim. names
    if(is.null(dim(beta_vals))) {
      # case 1: beta_vals is a vector 
      names(beta_vals) <- lamNames(lambda)
      } else {
        # case 2: beta_vals is a matrix
      colnames(beta_vals) <- lamNames(lambda)
    }
    
  }
  else beta_vals <- object$beta_vals[, which, drop = FALSE]
  
  if (drop){
    return(drop(beta_vals))
  } else{
    return(beta_vals)
  }
}
