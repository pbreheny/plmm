#' A function to implement singular value decomposition for a PLMM or LMM 
#' This is an internal function to \code{plmm_prep()}
#' 
#' @param X The standardized design matrix
#' @param k Optional integer argument indicating the number of singular values to use in a truncated SVD. See details. 
#' @param trunc Logical: should truncated SVD be used? 
#' @param trace Logical: should messages be printed to console? 
#' 
#' @details
#' The kind of SVD implemented here will depend on the combination of arguments supplied. 
#' (1) if trunc = FALSE
#'    use base::svd(K)
#' (2) if trunc = TRUE
#'    use RSpectra::svds(K, k = k)
#' 
#' @keywords internal

plmm_svd <- function(X, k, trunc, trace){
  # case 1: full SVD -----------------------------------  
  if(!trunc){
    if(trace){cat("\nUsing full SVD")}
    # NB: the is.null(S) keeps you from overwriting case 2 (a)
    decomp <- svd(X, nv = 0)
    d <- decomp$d
    U <- decomp$u
  } else {
    # case 2: truncated SVD -----------------------------
    if(trace){cat("\nUsing truncated SVD with k singular values")}
    decomp <- RSpectra::svds(A = X, nv = 0, k = k)
    d <- decomp$d
    U <- decomp$u
  }
  
  res <- list(d = d, U = U)
  return(res)
}