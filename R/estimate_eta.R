#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @export
#' @keywords internal
estimate_eta <- function(s, U, y, eta_star){
  # TODO: could also pass 'R' (the residuals... instead of y.)
  
  # coercion
  eta <- NULL
  
  # estimate eta 
  rot_y <- crossprod(U, y)
  opt <- stats::optimize(f=log_lik, c(0.01, 0.99), rot_y=rot_y, s=s)
  eta <- opt$minimum 
  
  
  return(eta)
}
