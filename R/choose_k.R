#' a function to choose k, the number of eigenvalues to use in truncated SVD
#' 
#' @param X can be either:
#'  *  the fully-imputed design matrix, or 
#'  *  a string specifying the path to an .rds object created by process_plink
#' @param start The starting number of eigenvalues. Defaults to floor(nrow(X)/10)
#' @param step The size (with respect to the number of observations) of the increments of increase in choosing k. Defaults to floor(nrow(X)/10).
#' @param eps The largest permissible difference between the true and approximate relatedness matrices. Defaults to 2.
#' @param trace Logical: should progress bars and messages be printed to the console? Defaults to TRUE. 
#' @param type The type of norm to use in determining distance between true and approximate K. Defaults to 'F', for Frobenious norm. See Matrix::norm() for details.
#' @param returnKapprox Logical: in addition to the list of SVD components for approximated K, should the approximation be returned as a matrix? Defaults to FALSE. 
#' @param returnK Logical: should the true K (as in relatedness_mat(X)) be returned? Defaults to FALSE.  
#' 
#' @return A list with at least 3 items:
#'  * K_svd: a list with the SVD components (d and u) for the approximated relatedness matrix K that can be passed to `plmm()`
#'  * k, the chosen number of eigenvalues 
#'  * delt, the distance between the true and approximated K matrices
#'  * optional: if returnKapprox, K_approx will be returned 
#'  * optional: if returnK, K will be returned 
#'  
#' @export
#'  
#' @examples
#' \dontrun{
#' choose_k(X = admix$X, start = 10, step = 10, trace = TRUE)
#' }
#'  
choose_k <- function(X,
                     start = NULL,
                     step = NULL,
                     eps = 2,
                     trace = TRUE,
                     type = "F",
                     returnKapprox = FALSE,
                     returnK = FALSE){
  
  # set defaults 
  if(is.null(start)){start <- floor(nrow(X)/10)}
  if(trace){cat("\nStarting k value is", start)}
  if(is.null(step)){step <- floor(nrow(X)/10)}
  if(trace){cat("\nStep size is", step)}
  
  # calculate true K 
  if(trace){cat('\nCalcuating the relatedness matrix')}
  std_X <- ncvreg::std(X)
  K <- penalizedLMM::relatedness_mat(std_X)
  
  # set up loop 
  k <- start
  it <- 1
  max_it <- (min(nrow(std_X), ncol(std_X)) - start)/step

  while(k < min(nrow(std_X), ncol(std_X)) | it < max_it){
    
    # truncated SVD 
    trunc <- RSpectra::svds(A = std_X, k = k, nv = 0)
    # calculate approximation
    S <- ((trunc$d)^2)/ncol(X)
    A_k <- trunc$u %*% tcrossprod(diag(S), trunc$u)
    
    # calculate difference 
    d <- abs(delta(K, A_k, type = type))
    
    if(d < eps){
      if(trace){cat("\nK approximation within specified epsilon reached at k =", k)}
      break
    }
    
    if(trace){cat('\nStepping up to next k value: k =', k+step)}
    k <- k + step
    it <- it + 1
  }
  
  # return list with approximated K 
  
  ret <- list(
    K_svd = trunc,
    k = k,
    delt = d
  )
  
  # optional items to return:
  if(returnK){ret$K_approx = A_k}
  if(returnKapprox){ret$K = K}
  
  return(ret)
  
}