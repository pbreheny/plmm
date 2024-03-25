#' A function to take the eigendecomposition of K
#' Note: This is faster than taking SVD of X when p >> n 
#'
#' @param std_X The *standardized* design matrix, stored as FBM.  
#' @param p The number of columns in the *unstandardized* design matrix
#' @param fbm_flag Logical: is std_X an FBM obejct? Passed from `plmm()`.
#'
#' @return A list with the eigenvectors and eigenvalues of K
#' @keywords internal
#'
eigen_K <- function(std_X, p, fbm_flag){
  # Note: std_X has already been scaled, so no need to do that here
  # calculate K (which will be stored in memory regardless of how std_X is stored)
  if(fbm_flag){
    K <- bigstatsr::big_tcrossprodSelf(std_X) 
    # TODO: this function returns an FBM. Need to think about adding the option 
    # for K to be stored filebacked. For now, make K stay in memory.
     K <- K[,]
  } else {
    K <- tcrossprod(std_X)
  }
  # take eigendecomposition
  decomp <- eigen(K)
  return(list(s = decomp$values/p,
              U = decomp$vectors,
              K = K))
}