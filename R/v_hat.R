#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix or list as returned by `choose_K()`
#' @return Vhat, a matrix representing the estimated variance 
#' 
#' @export
#' 
v_hat <- function(fit, K = NULL){
 
  # if K is supplied:
  if (!is.null(K)) {
    # case 1: K is from K_svd list 
    if(is.list(K)){
      K <- K$U %*% tcrossprod(diag(K$s), K$U)
      Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    } else if (is.matrix(K)){
    # case 2: K is a matrix 
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    }
  }
  
  # if K is not supplied:
  if(is.null(K) & is.list(fit)){
    K <- fit$U %*% tcrossprod(diag(fit$s), fit$U)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
  }
  
  
    return(Vhat)
}
 