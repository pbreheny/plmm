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
#' # colon, eg plmm:::plmm_prep()
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
                      trace = NULL, 
                      ...){
  
  
  ## coersion
  U <- s <- eta <- NULL
  
  # designate the dimensions of the original design matrix 
  n <- nrow(X)
  p <- ncol(X) 
  
  # standardize X
  # NB: the following line will eliminate singular columns (eg monomorphic SNPs)
  #  from the design matrix. 
  std_X <- ncvreg::std(X)
  
  # identify nonsingular values in the standardized X matrix  
  ns <- attr(std_X, "nonsingular")
  
  # keep only those penalty factors which penalize non-singular values 
  penalty.factor <- penalty.factor[ns]
  
  # designate the dimensions of the standardized design matrix, with only ns columns
  n_stdX <- nrow(std_X)
  p_stdX <- ncol(std_X)
  
  # name scaling and centering values (will need these later; see 'untransform()')
  std_X_details <- list(
    center = attr(std_X, 'center'), # singular columns have center = 0
    scale = attr(std_X, 'scale') # singular columns have scale = 0
  )
  
  # set default k and create indicator 'trunc' to pass to svd_X
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
  
  # handle the cases where no decomposition is needed: 
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
    s <- K$s # no need to adjust singular values by p; choose_k() does this via relatedness_mat()
    U <- K$U
  }

  # otherwise, need to do SVD:
  if(sum(c(flag1, flag2, flag3)) == 0){
    if(trace){cat("\nStarting decomposition.")}
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    # NB: relatedness_mat(X) standardizes X! 
    if(is.null(K) & is.null(s)){
      # NB: the is.null(s) keeps you from overwriting the 3 preceding special cases 
      
      # approach to decomposition:
      # n > p: take SVD of X
      # n <= p: construct K, then take eigen(K)
      if(n_stdX > p_stdX){
      if(trace){
        cat("\nCalculating the SVD of X")}
      svd_res <- svd_X(X = std_X, k = k, trunc = trunc, trace = trace)
      s <- (svd_res$d^2)*(1/p)
      U <- svd_res$U
      } else if (n_stdX <= p_stdX){
        if(trace){cat("\nCalculating eigendecomposition of K")}
        eigen_res <- eigen_K(std_X, p) 
        s <- eigen_res$s
        U <- eigen_res$U
      }
      
    } else {
      # last case: K is a user-supplied matrix 
      eigen_res <- eigen(K)
      s <- eigen_res$values
      U <- eigen_res$vectors
    }
    
  }

  # error check: what if the combination of args. supplied was none of the SVD cases above?
  if(is.null(s) | is.null(U)){
    stop("\nSomething is wrong in the SVD/eigendecomposition.
    \nThe combination of supplied arguments does not match any cases handled in 
         \n svd_X(), the internal function called by plmm() via plmm_prep().
         \n Re-examine the supplied arguments -- here are some common mistakes:
         \n \tDid you supply a list to K? Check its element names -- they must be 's' and 'U'. 
         \n \t \t *If you used choose_k(), make sure you are supplying the 'svd_K' element returned from that function's result as the 'K' here.*
         \n \tDid you intend to set diag_K = TRUE?
         \n \tDid you set diag_K = TRUE and specifiy a k value at the same time? This combination of arguments is incompatible.")
  }
  
  # estimate eta if needed
  if (is.null(eta_star)) {
    eta <- estimate_eta(s = s, U = U, y = y)$eta 
  } else {
    # otherwise, use the user-supplied value (this option is mainly for simulation)
    eta <- eta_star
  }
  
  
  # return values to be passed into plmm_fit(): 
  ret <- structure(list(
    p = p,
    n = n, 
    y = y,
    std_X = std_X,
    std_X_details = std_X_details,
    s = s,
    U = U,
    ns = ns,
    eta = eta, # carry eta over to fit 
    penalty.factor = penalty.factor,
    trace = trace,
    snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
