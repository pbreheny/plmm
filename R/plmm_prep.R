#' PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model
#' This is an internal function for \code{cv.plmm}
#'
#' @param std_X Column standardized design matrix. May include clinical covariates and other non-SNP data.
#' @param std_X_n The number of observations in std_X (integer)
#' @param std_X_p The number of features in std_X (integer)
#' @param p The number of features in the *original* design matrix X, including constant features
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized). This can also be a list, with components d and u (as returned by choose_k)
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `plmm()`. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Not used yet
#'
#' @return List with these components: 
#' * S: The singular values of K 
#' * U: the left singular values of K (same as left singular values of X). 
#' * snp_names: Formatted column names of the design matrix 
#'
#'@keywords internal
#'
#' @examples
#' 
#' \dontrun{
#' # this is an internal function; to call this, you would need to use the triple 
#' # colon, eg penalizedLMM:::plmm_prep()
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
#' prep2 <- plmm_prep(X = admix$X, y = admix$y, diag_K = TRUE, trace = TRUE)
#' }
#' 
plmm_prep <- function(std_X,
                      std_X_n,
                      std_X_p,
                      p,
                      y,
                      k = NULL,
                      K = NULL,
                      diag_K = NULL,
                      eta_star = NULL,
                      penalty.factor = rep(1, ncol(X)),
                      trace = NULL, 
                      ...){
  
  
  ## coersion
  U <- s <- eta <- NULL
  
  # set default k and create indicator 'trunc' to pass to plmm_svd
  if(is.null(k)){
    if('matrix' %in% class(std_X)){
      k <- min(std_X_n, std_X_p)
      trunc <- FALSE
    }
    
    if("FBM" %in% class(std_X)){
      # for FBM data, default k is 10% of the # of observations 
      k <- trunc(std_X_n*0.1)
      trunc <- TRUE
    }
    
  }
  if(!is.null(k) & (k < min(std_X_n, std_X_p))){
    trunc <- TRUE
  } else if(!(k %in% 1:min(std_X_n, std_X_p))){
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
    s <- rep(1, n)
    U <- diag(nrow = n)
  }
  # case 2: K is user-supplied diagonal matrix (like a weighted lm())
  flag2 <- diag_K & !is.null(K) & ('matrix' %in% class(K))
  if(flag2){
    if(trace){(cat("\nUsing supplied diagonal matrix for K, similar to a lm() with weights."))}
    s <- sort(diag(K), decreasing = T)
    U <- diag(nrow = n)[,order(diag(K), decreasing = T)]
  }
  # case 3: K is a user-supplied list, as passed from choose_k()
  flag3 <- !is.null(K) & ('list' %in% class(K))
  if(flag3){
    if(trace){cat("\nK is a list; will pass SVD components from list to model fitting.")}
    # case 3: K is a user-supplied list, as returned from choose_k()
    s <- K$s # no need to adjust singular values by p; choose_k() does this via relatedness_mat()
    U <- K$U
  }

  # otherwise, need to do SVD:
  if(trace){cat("\nStarting singular value decomposition.")}
  if(sum(c(flag1, flag2, flag3)) == 0){
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    if(is.null(K) & is.null(s)){
      # NB: the is.null(S) keeps you from overwriting the 3 preceding cases 
      if(trace){cat("\nUsing the default definition of the realized relatedness matrix.")}
      svd_res <- plmm_svd(X = std_X, k = k, trunc = trunc, trace = trace)
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
  ret <- structure(list(
    # include X and y here b/c I will need them for cross validation 
    std_X = std_X,
    y = y,
    s = s,
    U = U,
    trace = trace,
    snp_names = if (is.null(colnames(std_X))) paste("snp", 1:std_X_p, sep="") else std_X_p))
  
  return(ret)
  
  
  
}
