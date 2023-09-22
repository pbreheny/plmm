#' PLMM prep FBM: a function to run checks, SVD, and rotation prior to fitting a PLMM model when X is an FBM
#'
#' @param X Design matrix of type Filebacked Big Matrix (FBM). May include clinical covariates and other non-SNP data.
#' @param meta A list with the appropriate meta-data to accompany X. This will include 'map' and 'fam' data frames, as well as an 'ns' vector. See `process_plink()` and `get_data()` for details.
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized). This can also be a list, with components d and u (as returned by choose_k)
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `plmm()`. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Other arguments to `bigstatsr::big_randomSVD()`
#'
#' @return List with these components: 
#' * ncol_X: The number of columns in the original design matrix 
#' * std_X: standardized design matrix 
#' * y: The vector of outcomes 
#' * S: The singular values of K 
#' * U: the left singular values of K (same as left singular values of X). 
#' * ns: the indices for the nonsingular values of std_X
#' * penalty.factor: the penalty factors for the penalized non-singular values 
#' * snp_names: Formatted column names of the design matrix 
#'
#'@keywords internal
plmm_prep_fbm <- function(X,
                          meta,
                          y,
                          k = NULL,
                          K = NULL,
                          diag_K = NULL,
                          eta_star = NULL,
                          penalty.factor = rep(1, ncol(X)),
                          trace = FALSE, 
                          ...){
  
  ## coersion
  U <- S <- SUX <- SUy <- eta <- NULL
  
  # designate the dimensions of the design matrix 
  p <- X$nrow 
  n <- X$ncol
  
  # note what values are non-singular
  ns <- meta$ns
  
  # set default k 
  if(is.null(k)){
    k <- min(n,p)
  }
  
  if(is.null(diag_K)){
    diag_K <- FALSE
  } else {
    diag_K <- TRUE
  }
  
  # keep only those penalty factors which penalize non-singular values 
  penalty.factor <- penalty.factor[ns]

  # designate the dimensions of the standardized design matrix, with only ns columns
  n_stdX <- X$nrow
  p_stdX <- length(ns)
  
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
    S <- K$d # no need to adjust singular values by p; choose_k() does this via relatedness_mat()
    U <- K$u
  }
  
  # otherwise, need to do SVD:
  if(trace){cat("\nStarting singular value decomposition.")}
  if(sum(c(flag1, flag2, flag3)) == 0){
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    # NB: relatedness_mat(X) standardizes X! 
    if(is.null(K) & is.null(S)){
      # NB: the is.null(S) keeps you from overwriting the 3 preceding cases 
      if(trace){cat("\nUsing the default definition of the realized relatedness matrix.")}
      K <- relatedness_mat_fbm(X, ns)
    }
    
    svd_res <- bigstatsr::big_randomSVD(X = K,
                                        k = k,
                                        ...)
    S <- svd_res$d
    U <- svd_res$u |> as_FBM()
  }
  
  # error check: what if the combination of args. supplied was none of the SVD cases above?
  if(is.null(S) | is.null(U)){
    stop("\nSomething is wrong in the SVD. The combination of supplied arguments does not match any cases handled in 
         \n plmm_svd(), the internal function called by plmm() via plmm_prep().
         \n Re-examine the supplied arguments.")
  }
  
  # return values to be passed into plmm_fit(): 
  ret <- structure(list(ncol_X = n,
                        nrow_X = p, 
                        X = X,
                        y = y,
                        # NB: watch the line below; different from plmm_prep()
                        center = ifelse(exists('decomp'), decomp$center, meta$center),
                        scale = ifelse(exists('decomp'), decomp$scale, meta$scale),
                        S = S,
                        U = U,
                        ns = meta$ns,
                        eta = eta_star, # carry eta over to fit 
                        penalty.factor = penalty.factor,
                        trace = trace,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
