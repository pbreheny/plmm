#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param S The diagonal matrix whose nonzero values are the singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @importFrom zeallot %<-%
#' @export
#' 
estimate_eta <- function(S, U, y, eta_star){
  # TODO: could also pass 'R' (the residuals... instead of y.)
  
  # coercion
  eta <- NULL
  
  # estimate eta 
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=logLik, c(0.01, 0.99), Uy=Uy, S=S)
  eta <- opt$minimum 
  
  
  return(eta)
}
