#' A function to implement singular value decomposition for a PLMM or LMM 
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
#' In case 1 (K is NULL): 
#'  * (a) if diag_K = TRUE: K is treated as the identity matrix. 
#'     All S = rep(1, n), and U = diag(n)
#'  * (b) diag_K = FALSE and k = min(n,p): 
#'    use base::svd(K), where K = relatedness_mat(std_X)
#'  * (c) diag_K = FALSE and k < min(n,p): 
#'    use RSpectra::svds(K, k = k), where K = relatedness_mat(std_X)
#'  
#' In case 2 (K is a matrix): 
#'  * (a): if diag_K = TRUE: observations will be treated as unrelated, 
#'    where S = sort(diag(K), decreasing = T) and columns of U are the columns of an 
#'    identity matrix sorted on the values of S.  
#'  * (b): if k = min(n,p): use base::svd(K)
#'  * (c): if k < min(n,p): use RSpectra::svds(K, k)
#'  
#' In case 3: 
#' if K is a list, D and U are simply passed from \code{choose_k()} & D is transformed to S. 
#' 
#' Any scenarios outside of those outlined will error out. 
#' 
#' @keywords internal

plmm_svd <- function(std_X, n, p, diag_K, K, k, trace){

  # prelim step: due to the presence of singular values, we need to distinguish 
  #   between the dimensions of X and the dimensions of std_X when choosing k: 
  n_stdX <- nrow(std_X)
  p_stdX <- ncol(std_X)

  # case 1: K is not specified ------------------------------
  if(is.null(K)){
    ## case 1 (a): if K is the identity matrix, then no need for SVD!
    if(diag_K){
      if(trace){(cat("Using identity matrix for K."))}
      S <- rep(1, n)
      U <- diag(nrow = n)
    } else if (!diag_K){
      ## case 1 (b): K is not specified, diag_K = FALSE, and k = min(n,p)
      if(trace){cat("No K specified - will use default definition of the \n realized relatedness matrix.\n")}
      K <- relatedness_mat(std_X) 
      # TODO: thinking about how to write relatedness_mat() in C...
      
      # remember: if I want all the singular values (which is k = min(n,p)), use base::svd
      if(k == min(n_stdX,p_stdX)){
        if(trace){cat("Using all singular values (full SVD)")}
        decomp <- svd(K, nv = 0)
      }
      ## case 1 (c): K is not specified, not diagonal, and k < min(n,p)
      # if I want fewer singular values than min(n,p), use RSpectra decomposition method:
      if(trace){cat('Using truncated SVD with # of singular values specified by k')}
      if (k < min(n_stdX,p_stdX)){
        decomp <- RSpectra::svds(A = K, nv = 0, k = k)
      }
      
      S <- decomp$d # relatedness_mat() internally scales by 1/p, so no need for that here
      U <- decomp$u 
      
    }
    # case 2: K is a user-specified matrix -----------------------------------  
  } else if (!is.null(K) & 'matrix' %in% class(K)){
    
    S <- U <- NULL
    
    ## case 2 (a): K is a diagonal matrix and diag_K = TRUE
    if(diag_K){
      if(trace){(cat("Using diagonal matrix for K, equivalent to a lm() with weights."))}
      S <- sort(diag(K), decreasing = T)
      U <- diag(nrow = n)[,order(diag(K), decreasing = T)]
    }
    
    # case 2 (b): K is a matrix, but not diagonal, and k = min(n,p)
    if(k == min(n_stdX, p_stdX) & is.null(S)){
      # NB: the is.null(S) keeps you from overwriting case 2 (a)
      decomp <- svd(K, nv = 0)
      S <- decomp$d/p # need to scale singular values by 1/p here!
      U <- decomp$u
    } 
    # case 2 (c): K is a matrix and k < min(n,p)
    if(k < min(n_stdX,p_stdX) & is.null(S)){
      decomp <- RSpectra::svds(A = K, nv = 0, k = k)
      S <- decomp$d/p # need to scale singular values by 1/p here!
      U <- decomp$u
    }
    # case 3: K is a list ---------------------------------------------------
  } else if(!is.null(K) & is.list(K)){
    # case 3: K is a user-supplied list, as returned from choose_k()
    S <- K$d # no need to adjust singular values by p; choose_k() does this via relatedness_mat()
    U <- K$u
    
  }
  
  # error check: what if the combination supplied was none of the above?
  if(is.null(S) | is.null(U)){
    stop("Something is wrong in the SVD. The combination of supplied arguments does not match any cases handled in 
         \n plmm_svd(), the internal function called by plmm() via plmm_prep().
         \n Re-examine the supplied arguments -- should you have set diag_K = TRUE?
         \n Or did you set diag_K = TRUE and specifiy a k value at the same time?")
  }
  
  res <- list(S = S, U = U)
  return(res)
}