#' PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model
#' This is an internal function for \code{cv.plmm}
#'
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized). This can also be a list, with components d and u (as returned by choose_k)
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `plmm()`. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Not used yet
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
plmm_prep <- function(X,
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
  
  # keep only those penalty factors which penalize non-singular values 
  penalty.factor <- penalty.factor[ns]
  
  # designate the dimensions of the design matrix 
  p <- ncol(X) 
  n <- nrow(X)
  # calculate SVD
  if(!(k %in% 1:min(n,p))){stop("k value is out of bounds. \nIf specified, k must be in the range from 1 to min(nrow(X), ncol(X))")}
  ## case 1: K is not specified (default to realized relatedness matrix)
  # check: if K is diagonal, then no need for SVD! 
  if(is.null(K) & diag_K){
    if(trace){(cat("Using diagonal for K, so observations are treated as unrelated."))}
    S <- rep(1, n)
    U <- diag(n)
  } else if (is.null(K) & !diag_K){
    if(trace){cat("No K specified - will use default definition of the \n realized relatedness matrix.\n")}
    
    # if I want all the singular values (which is k = min(n,p)), use base::svd
    if(k == min(nrow(std_X),ncol(std_X))){
      decomp <- svd(std_X, nv = 0)
    }
    # otherwise, if I want fewer singular values than min(n,p), use RSpectra decomposition method:
    if (k < min(nrow(std_X),ncol(std_X))){
      decomp <- RSpectra::svds(A = std_X, nv = 0, k = k)
    }
    
    D <- decomp$d
    U <- decomp$u
    S <- (D^2)/p # singular values of K, the realized relationship matrix
    
  } else if (!is.null(K) & 'matrix' %in% class(K)){
    ## case 2: K is a user-specified matrix
    S <- U <- NULL
    
    # again, decomposition depends on choice of k
    if(k == min(nrow(std_X),ncol(std_X))){
      decomp <- svd(K, nv = 0)
    } 
    
    if(k < min(nrow(std_X),ncol(std_X))){
      decomp <- RSpectra::svds(A = K, nv = 0, k = k)
    }
    S <- decomp$d # if K matrix was user-specified, then no need to transform D here
    U <- decomp$u
  } else if(!is.null(K) & is.list(K)){
    # case 3: K is a user-supplied list, as returned from choose_k()
      S <- ((K$d)^2)/p
      U <- K$u

  }

  # return values to be passed into plmm_fit(): 
  ret <- structure(list(ncol_X = ncol(X),
                        nrow_X = nrow(X), 
                        y = y,
                        std_X = std_X,
                        S = S,
                        U = U,
                        ns = ns,
                        eta = eta_star, # carry eta over to fit 
                        penalty.factor = penalty.factor,
                        trace = trace,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
