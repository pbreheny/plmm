#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix or list as returned by `choose_K()`
#' @return Vhat, a matrix representing the estimated variance 
#' 
#' @export
#' 
v_hat <- function(fit, K = NULL){
  
  if(is.null(fit) & is.null(K))stop("Either fit or K must be specified")
  # if K is supplied:
  if (!is.null(K)) {
    # case 1: K is from K_svd list 
    if(is.list(K)){
      # TODO: see if I can find a better computational approach here -- 
      U <- K$U
      s <- K$s
      K <- K$U %*% tcrossprod(diag(s), K$U)
      Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    } else if (is.matrix(K)){
    # case 2: K is a matrix 
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    }
  }

  # if K is not supplied:
  if(is.null(K) & is.list(fit)){
    # case 3: K is default (RRM)
    U <- fit$U
    s <- fit$s
    K <- U %*% tcrossprod(diag(s), U)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
  }
  
  
    return(Vhat)
}
 