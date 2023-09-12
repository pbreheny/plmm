#' a function to prepare data for an *unpenalized* LMM 
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized).
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `lmm()`. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @param ... Not used yet
#'
#' @return List with these components: 
#' * ncol_X: The number of columns in the original design matrix 
#' * std_X: standardized design matrix 
#' * y: The vector of outcomes 
#' * S: The singular values of K 
#' * U: the left singular values of K (same as left singular values of X). 
#' * ns: the indices for the nonsingular values of std_X
#' * snp_names: Formatted column names of the design matrix 
#'
#'@keywords internal
lmm_prep <- function(X,
                      y,
                      k = NULL,
                      K = NULL,
                      diag_K = NULL,
                      eta_star = NULL,
                      returnX = TRUE,
                      trace = FALSE, ...){
  
  
  ##coercion
  U <- S <- SUX <- SUy <- eta <- NULL
  
  # standardize X
  # NB: the following line will eliminate singular columns (eg monomorphic SNPs)
  #  from the design matrix. 
  std_X <- ncvreg::std(X) 
  
  # set default k 
  if(is.null(k)){
    k <- min(nrow(std_X),ncol(std_X))
  }
  
  if(is.null(diag_K)){
    diag_K <- FALSE
  } else {
    diag_K <- TRUE
  }
  
  # identify nonsingular values in the standardized X matrix  
  ns <- attr(std_X, "nonsingular")
  
  # designate the dimensions of the design matrix 
  p <- ncol(X) 
  n <- nrow(X)
  
  # calculate SVD
  if(!(k %in% 1:min(n,p))){stop("k value is out of bounds. \nIf specified, k must be in the range from 1 to min(nrow(X), ncol(X))")}
  svd_res <- plmm_svd(std_X, n, p, diag_K, K, k, trace)
  
  # return values to be passed into lmm_fit(): 
  ret <- structure(list(ncol_X = ncol(X),
                        nrow_X = nrow(X), 
                        y = y,
                        std_X = std_X,
                        S = svd_res$S,
                        U = svd_res$U,
                        ns = ns,
                        eta = eta_star, # carry eta over to fit 
                        trace = trace,
                        returnX = returnX,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
