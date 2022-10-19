#' Estimate eta using a null Gaussian penalizedLMM model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param X The matrix of SNP data. If K is not supplied, X is required. 
#' @param y Continuous outcome vector. If not supplied, K is treated as known and eta is not estimated.
#' @param K Estimated or known similarity matrix. By default, K is the realized relationship matrix, \eqn{\frac{1}{p}XX^T}, where \eqn{p} is the number of columns in X
#' @param k An integer between 1 and \code{nrow(K)} indicating the number of singular values requested *if* package \code{RSpectra} is installed. Defaults to NULL. 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples
#' res <- plmm_null(y = admix$y, X = admix$X)
#' (res$eta) # estimate of eta 

plmm_null <- function(y, X = NULL, K = NULL, k = NULL){
  
  # error checking 
  if(is.null(K) & is.null(X)){stop("If K is not supplied, X is required")}

  # if K is null, calculate realized relationship matrix
  if(is.null(K) & !is.null(X)){
    K <- relatedness_mat(X)
  } 
  
  # estimate eta 
  S <- U <- NULL
  
  if("RSpectra" %in% rownames(installed.packages()) & !is.null(k)){
    svd_K <- RSpectra::svds(K, nv = 0, k = k)
    # NB: nv=0 avoids the calculation of right singular vectors, which we will not use
  } else {
    svd_K <- svd(K)
  }
  
  S <- svd_K$d
  U <- svd_K$u
 
  
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=logLik, c(0.01, 0.99), Uy=Uy, S=S)
  eta <- opt$minimum
  return(list("S" = S, "U" = U, "eta" = eta))
}
