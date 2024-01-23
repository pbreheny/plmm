#' A function to implement singular value decomposition for a PLMM or LMM 
#' This is an internal function to \code{plmm_prep()}
#' 
#' @param std_X The standardized design matrix, which could be stored as an FBM object or in memory. 
#' @param k Optional integer argument indicating the number of singular values to use in a truncated SVD. See details. 
#' @param trunc Logical: should truncated SVD be used? 
#' @param fbm_flag: Logical: is std_X and FBM? Passed from `plmm()`.
#' @param trace Logical: should messages be printed to console? 
#' @param ... Additional arguments to bigstatsr::big_randomSVD()
#' @details
#' The kind of SVD implemented here will depend on the combination of arguments supplied. 
#' (1) if trunc = FALSE
#'    use base::svd(K)
#' (2) if trunc = TRUE
#'    use RSpectra::svds(K, k = k)
#' 
#' @keywords internal

svd_X <- function(std_X, k, trunc, fbm_flag, trace, ...){
  if(fbm_flag){
    if(!trunc){
      # case 1: full SVD ----------------------------------- 
      if(trace){cat("\nUsing full SVD, which requires loading std_X into memory.")}
      decomp <- svd(std_X[,], nv = 0)
      d <- decomp$d
      U <- decomp$u
    } else {
      # case 2: truncated SVD -----------------------------
      if(trace){cat("\nUsing truncated SVD with", k ,"singular values")}
      decomp <- bigstatsr::big_randomSVD(X = X, k = k, ...)
      d <- decomp$d
      U <- decomp$u |> bigstatsr::as_FBM()
    }
  } else {
    # case 1: full SVD -----------------------------------  
    if(!trunc){
      if(trace){cat("\nUsing full SVD")}
      decomp <- svd(X, nv = 0)
      d <- decomp$d
      U <- decomp$u
    } else {
      # case 2: truncated SVD -----------------------------
      if(trace){cat("\nUsing truncated SVD with", k ,"singular values")}
      decomp <- RSpectra::svds(A = X, k = k, ...)
      d <- decomp$d
      U <- decomp$u 
    }
  }
  
  res <- list(d = d, U = U)
  return(res)
}