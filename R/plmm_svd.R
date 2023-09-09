#' A function to implement singular value decomposition for a PLMM 
#' This is an internal function to \code{plmm_prep()}
#' 
#' @param std_X The standardized X matrix
#' @param n The number of rows in X
#' @param p The number of columns in X 
#' @param diag_K Logical: is K a diagonal matrix? Defaults to FALSE. 
#' @param K Optional argument indicating the relatedness matrix to be used. See details. 
#' @param k Optional integer argument indicating the number of singular values to use in a truncated SVD. See details. 
#' @param trace Logical: should messages be printed to console? 
#' 
#' @details
#' The kind of SVD implemented here will depend on the combination of arguments supplied. 
#' Let svd(K) be defined as t(U)SU, since K is symmetric in a PLMM context. We can
#' let svd(std_X) be defined as t(U)DU, and write S = (D^2)/p (a good derivation exercise, for those who are interested). 
#' 
#' Broadly speaking, there are three cases: 
#' (1) K is not user-specified, (2) K is a matrix, and (3) K is a list. 
#' 
#' In case 1: 
#'  * (a) if K is not specified but diag_K = TRUE: observations will be treated as unrelated, 
#'    where S = sort(diag(K), decreasing = T) and columns of U are the columns of an 
#'    identity matrix sorted on the values of S.  
#'  * (b) if K is not specified, diag_K = FALSE, and k = min(n,p): 
#'    use base::svd(K), where K = relatedness_mat(std_X)
#'  * (c) if K is not specified, diag_K = FALSE, and k < min(n,p): 
#'    use RSpectra::svds(K, k = k), where K = relatedness_mat(std_X)
#'  
#' In case 2: 
#'  * (a): if K is a matrix and k = min(n,p): use base::svd(K)
#'  * (b): if K is a matrix and k < min(n,p): use RSpectra::svds(K, k)
#'  
#' In case 3: 
#' if K is a list, D and U are simply passed from \code{choose_k()} & D is transformed to S. 
#' 
#' Any scenarios outside of these 3 will error out. 
#' 
#' @keywords internal

plmm_svd <- function(std_X, n, p, diag_K, K, k, trace){
  
  # prelim step: due to the presence of singular values, we need to distinguish 
  #   between the dimensions of X and the dimensions of std_X when choosing k: 
  n_stdX <- nrow(std_X)
  p_stdX <- ncol(std_X)

  # case 1: K is not specified (default to realized relatedness matrix)
  if(is.null(K) & diag_K){
    ## case 1 (a): if K is diagonal, then no need for SVD! 
    if(trace){(cat("Using diagonal for K, so observations are treated as unrelated."))}
    S <- sort(diag(K), decreasing = T)
    U <- diag(nrow = n)[,order(diag(K), decreasing = T)]
  } else if (is.null(K) & !diag_K){
    ## case 1 (b): K is not specified, not diagonal, and k = min(n,p)
    if(trace){cat("No K specified - will use default definition of the \n realized relatedness matrix.\n")}
    K <- relatedness_mat(std_X) # TODO: write this function in C 
    # remember: if I want all the singular values (which is k = min(n,p)), use base::svd
    if(k == min(n_stdX,p_stdX)){
      decomp <- svd(K, nv = 0)
    }
    ## case 1 (c): K is not specified, not diagonal, and k < min(n,p)
    # if I want fewer singular values than min(n,p), use RSpectra decomposition method:
    if (k < min(n_stdX,p_stdX)){
      decomp <- RSpectra::svds(A = K, nv = 0, k = k)
    }
    
    S <- decomp$d
    U <- decomp$u # NB: we need this -1 multiplier!! 
    
  } else if (!is.null(K) & 'matrix' %in% class(K)){
    # case 2: K is a user-specified matrix
    S <- U <- NULL
    
    # case 2 (a): K is a matrix and k = min(n,p)
    if(k == min(n_stdX, p_stdX)){
      decomp <- svd(K, nv = 0)
    } 
    # case 2 (b): K is a matrix and k < min(n,p)
    if(k < min(n_stdX,p_stdX)){
      decomp <- RSpectra::svds(A = K, nv = 0, k = k)
    }
    S <- decomp$d # if K matrix was user-specified, then no need to transform D here
    U <- decomp$u
  } else if(!is.null(K) & is.list(K)){
    # case 3: K is a user-supplied list, as returned from choose_k()
    S <- K$d
    U <- K$u
    
  }
  res <- list(S = S, U = U)
  return(res)
}