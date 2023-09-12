#' Coef method for "lmm" class
#'
#' @param object An object of class "lmm."
#' @param which Vector of indices for which coefficients to return.
#' @param ... Additional arguments.
#' 
#' @rdname coef.lmm
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' fit <- lmm(X = pedigree$X, y = pedigree$clinical$y, K = pedigree$K)
#' coef.lmm(fit)
#' }
#' 


coef.lmm <- function(object, which = 1:length(object$beta_vals), ...){
  beta_vals <- object$beta_vals[which]
  return(beta_vals)
  
}

# TODO: add summary & print.summary methods here