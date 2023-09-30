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
  
  
  
  ## coersion
  U <- S <- SUX <- SUy <- eta <- NULL
  
  # designate the dimensions of the original design matrix 
  n <- nrow(X)
  p <- ncol(X) 
  
  # standardize X
  # NB: the following line will eliminate singular columns (eg monomorphic SNPs)
  #  from the design matrix. 
  std_X <- ncvreg::std(X)
  
  # identify nonsingular values in the standardized X matrix  
  ns <- attr(std_X, "nonsingular")
  
  # designate the dimensions of the standardized design matrix, with only ns columns
  n_stdX <- nrow(std_X)
  p_stdX <- ncol(std_X)
  
  # set default k and create indicator 'trunc' to pass to plmm_svd
  if(is.null(k)){
    k <- min(n_stdX, p_stdX)
    trunc <- FALSE
  } else if(!is.null(k) & (k < min(n_stdX, p_stdX))){
    trunc <- TRUE
  } else if(!(k %in% 1:min(n_stdX,p_stdX))){
    stop("\nk value is out of bounds. 
         \nIf specified, k must be an integer in the range from 1 to min(nrow(X), ncol(X)). 
         \nwhere X does not include singular columns. For help detecting singularity,
         \nsee ncvreg::std()")
  }
  
  # set default: if diag_K not specified, set to false
  if(is.null(diag_K)){diag_K <- FALSE}
  
  # handle the cases where no SVD is needed: 
  # case 1: K is the identity matrix 
  flag1 <- diag_K & is.null(K)
  if(flag1){
    if(trace){(cat("\nUsing identity matrix for K."))}
    S <- rep(1, n)
    U <- diag(nrow = n)
  }
  # case 2: K is user-supplied diagonal matrix (like a weighted lm())
  flag2 <- diag_K & !is.null(K) & ('matrix' %in% class(K))
  if(flag2){
    if(trace){(cat("\nUsing supplied diagonal matrix for K, similar to a lm() with weights."))}
    S <- sort(diag(K), decreasing = T)
    U <- diag(nrow = n)[,order(diag(K), decreasing = T)]
  }
  # case 3: K is a user-supplied list, as passed from choose_k()
  flag3 <- !is.null(K) & ('list' %in% class(K))
  if(flag3){
    if(trace){cat("\nK is a list; will pass SVD components from list to model fitting.")}
    # case 3: K is a user-supplied list, as returned from choose_k()
    S <- K$s # no need to adjust singular values by p; choose_k() does this via relatedness_mat()
    U <- K$U
  }
  
  # otherwise, need to do SVD:
  if(trace){cat("\nStarting singular value decomposition.")}
  if(sum(c(flag1, flag2, flag3)) == 0){
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    # NB: relatedness_mat(X) standardizes X! 
    if(is.null(K) & is.null(S)){
      # NB: the is.null(S) keeps you from overwriting the 3 preceding cases 
      if(trace){cat("\nUsing the default definition of the realized relatedness matrix.")}
      svd_res <- plmm_svd(X = X, k = k, trunc = trunc, trace = trace)
      s <- (svd_res$d^2)*(1/p)
      U <- svd_res$U
    } else {
      # last case: K is a user-supplied matrix
      svd_res <- plmm_svd(X = K, k = k, trunc = trunc, trace = trace)
      s <- svd_res$d
      U <- svd_res$U
    }
    
  }
  
  # error check: what if the combination of args. supplied was none of the SVD cases above?
  if(is.null(s) | is.null(U)){
    stop("\nSomething is wrong in the SVD. The combination of supplied arguments does not match any cases handled in 
         \n plmm_svd(), the internal function called by plmm() via plmm_prep().
         \n Re-examine the supplied arguments -- should you have set diag_K = TRUE?
         \n Or did you set diag_K = TRUE and specifiy a k value at the same time?")
  }
  
  
  
  # return values to be passed into plmm_fit(): 
  ret <- structure(list(p= p,
                        n = n, 
                        y = y,
                        std_X = std_X,
                        S = S,
                        U = U,
                        ns = ns,
                        eta = eta_star, # carry eta over to fit 
                        trace = trace,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
