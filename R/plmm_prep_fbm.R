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
  # TODO: Think about implementing standardization here with big_apply()
  # std_X <- bigstatsr::big_apply(X = X,
  #                               a.FUN = std_fbm,
  #                               a.combine = cbind,
  #                               ncores = bigstatsr::nb_cores())
  # browser()
  
  # designate the dimensions of the design matrix 
  p <- meta$X$nrow 
  n <- meta$X$ncol
  
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
  # browser()
  # calculate SVD
  if(!(k %in% 1:min(n,p))){stop("k value is out of bounds. \nIf specified, k must be in the range from 1 to min(nrow(X), ncol(X))")}
  ## case 1: K is not specified (default to realized relatedness matrix)
  # check: if K is diagonal, then no need for SVD! 
  if(is.null(K) & diag_K){
    if(trace){(cat("Using diagonal for K, so observations are treated as unrelated."))}
    S <- rep(1, n)
    U <- diag(n)
  } else 
    # TODO: double check, but I don't think the scenario below will apply when X is an FBM
    # if (is.null(K) & !diag_K){
    # if(trace){cat("No K specified - will use default definition of the \n realized relatedness matrix.\n")}
    # 
    # # if I want all the singular values (which is k = min(n,p)), use base::svd
    # if(k == min(nrow(std_X),ncol(std_X))){
    #   decomp <- svd(std_X, nv = 0)
    # }
    # otherwise, if I want fewer singular values than min(n,p), use RSpectra decomposition method:
    if (k < min(n,p)){
      decomp <- bigstatsr::big_randomSVD(X,
                              # standardization happens because of the line below
                              fun.scaling = bigstatsr::as_scaling_fun(center.col = meta$center, scale.col = meta$scale),
                              k = k,
                              ind.col = ns,
                              ...)

    
    D <- decomp$d
    U <- decomp$u
    S <- (D^2)/p # singular values of K, the realized relationship matrix
    
  } else if (!is.null(K) & 'matrix' %in% class(K)){
    ## case 2: K is a user-specified matrix
    S <- U <- NULL
    # 
    # # again, decomposition depends on choice of k
    # if(k == min(n,p)){
    #   decomp <- svd(K, nv = 0)
    # } 
    # 
    if(k < min(n,p)){
      decomp <- bigstatsr::big_randomSVD(K,
                              fun.scaling = bigstatsr::big_scale(center = TRUE, scale = TRUE),
                              k = k,
                              ...)
    }
    S <- decomp$d # if K matrix was user-specified, then no need to transform D here
    U <- decomp$u
  } else if(!is.null(K) & is.list(K)){
    # case 3: K is a user-supplied list, as returned from choose_k()
    S <- ((K$d)^2)/p
    U <- K$u
    
  }
  
  # return values to be passed into plmm_fit(): 
  ret <- structure(list(ncol_X = n,
                        nrow_X = p, 
                        X = X,
                        y = y,
                        # NB: watch the line below; different from plmm_prep()
                        std_X = ifelse(exists('decomp'),
                                       list(center = decomp$center,
                                            scale = decomp$scale),
                                       list(center = meta$center,
                                            scale = meta$scale)),
                        S = S,
                        U = U,
                        ns = meta$ns,
                        eta = eta_star, # carry eta over to fit 
                        penalty.factor = penalty.factor,
                        trace = trace,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
