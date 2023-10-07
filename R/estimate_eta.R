#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @export
#' @keywords internal
estimate_eta <- function(s, U, y){
  
  # coercion
  eta <- NULL
  # estimate eta 
  if("FBM" %in% class(U)){
    n <- U$nrow
    rot_y <- bigstatsr::big_cprodVec(U, y)
    opt <- stats::optimize(f=log_lik,
                           interval = c(0.01, 0.99),
                           rot_y=rot_y,
                           s=s,
                           n=n)
    eta <- opt$minimum 
    
  } else {
    n <- nrow(U)
    rot_y <- crossprod(U, y)
    opt <- stats::optimize(f=log_lik,
                           interval = c(0.01, 0.99),
                           rot_y=rot_y,
                           s=s,
                           n=n)
    eta <- opt$minimum 
    
  }
  
  
  return(eta)
}
