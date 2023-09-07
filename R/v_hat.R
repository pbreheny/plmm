#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional list or matrix as returned by `choose_K()`
#' @return Vhat, a matrix representing the estimated variance 
#' 
#' @export
#' 
v_hat <- function(fit, K = NULL){
  
  if(is.null(fit) & is.null(K))stop("Either fit or K must be specified")
  
  # if K is supplied:
  if (!is.null(K) & is.list(K)) {
    # case 1: K is from K_svd list 
    K <- K$U %*% tcrossprod(diag(K$S), K$U)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$nrow_X)
  } else if (!is.null(K) & is.matrix(K)){
    # case 2: K is a matrix (either from returnK or returnKapprox)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$nrow_X)
  }
  
  # if K is not supplied:
  if(is.null(K) & is.list(fit)){
    K <- fit$U %*% tcrossprod(diag(fit$S), fit$U)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$nrow_X)
  }
  
  
    return(Vhat)
}
 