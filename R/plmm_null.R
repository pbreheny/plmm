#' Estimate eta using a null Gaussian penalizedLMM model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param X The matrix of SNP data. If K is not supplied, X is required. 
#' @param y Continuous outcome vector. If not supplied, V is treated as known and eta is not estimated.
#' @param K Estimated or known similarity matrix. By default, K is the realized relationship matrix, \eqn{\frac{1}{p}XX^T}, where \eqn{p} is the number of columns in X
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples
#' res <- plmm_null(y = admix$y, X = admix$X)
#' (res$eta) # estimate of eta 

plmm_null <- function(y, X = NULL, K = NULL){
  
  # error checking 
  if(is.null(K) & is.null(X)){stop("If K is not supplied, X is required")}

  # if K is null, calculate realized relationship matrix
  if(is.null(K) & !is.null(X)){
    p <- ncol(X)
    K <- tcrossprod(X, X) * (1/p)
  } 
  
  # estimate eta 
  S <- U <- NULL
  c(S, U) %<-% svd(K)[1:2]
  # NB: svd() returns components d, u, and v (in that order)
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=plmm_nll, c(0.01, 0.99), Uy=Uy, S=S)
  eta <- opt$minimum
  return(list(S=S, U=U, eta=eta))
}
