#' a function to choose k, the number of eigenvalues to use in truncated SVD
#' 
#' @param X can be either:
#'  *  the fully-imputed design matrix, or 
#'  *  a string specifying the path to an .rds object created by process_plink
#' @param start The starting number of eigenvalues. Defaults to floor(nrow(X)/3)
#' @param step The size (with respect to the number of observations) of the increments of increase in choosing k
#' @param eps The largest permissible difference between the true and approximate relatedness matrices. Defaults to 0.01.
#' @param trace Logical: should progress bars and messages be printed to the console? Defaults to TRUE. 
#' @return A list with two items:
#'  * An approximated relatedness matrix K that can be passed to `plmm()`
#'  * k, the chosen number of eigenvalues 


#' a helper function to calcuate the 'distance' between 2 matrices 
delta <- function(A, B, type = "F"){
  Matrix::norm(A, type = type) - Matrix::norm(B, type = type)
}

choose_k <- function(X, start = NULL, step = NULL, eps = 0.01, trace = TRUE){

  # set defaults 
  if(is.null(start)){start <- floor(nrow(X)/3)}
  if(trace){cat("Starting k value is ", start)}
  if(is.null(step)){step <- floor(nrow(X))/10}
  if(trace){cat("Step size is ", step)}
  
  # calculate true K 
  K <- penalizedLMM::relatedness_mat(X)

  # set up loop 
  k <- start
  if(trace){pb <- txtProgressBar(min = 0, max = k, style = 3)}
  while(k < min(nrow(X), ncol(X))){
    
    # calculate approximation
    trunc <- RSpectra::svds(A = K, k = k, nv = 0)
    A_k <- trunc$u %*% tcrossprod(diag(trunc$d), trunc$u)
    
    # calculate difference 
    d <- abs(delta(K, A_k))
    
    if(d < eps){break}
    
    k <- k + step
    if(trace){setTxtProgressBar(pb, k)}
    
  }
  
  # return approximated K 
  if(trace){cat("K approximation within specified epsilon reached at k = ", k)}
  return(list(K_approx = A_k,
              k = k))

}
