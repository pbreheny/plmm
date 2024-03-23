#' PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model
#' This is an internal function for \code{cv.plmm}
#'
#' @param std_X Column standardized design matrix. May include clinical covariates and other non-SNP data.
#' @param std_X_n The number of observations in std_X (integer)
#' @param std_X_p The number of features in std_X (integer)
#' @param n The number of instances in the *original* design matrix X. This should not be altered by standardization.
#' @param p The number of features in the *original* design matrix X, including constant features
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized). This can also be a list, with components d and u (as returned by choose_k)
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `plmm()`. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param fbm_flag Logical: is std_X an FBM type object? This is set internally by `plmm()`. 
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Not used yet
#'
#' @returns List with these components: 
#' * n: the number of rows in the original design matrix
#' * p: the number of columns in the original design matrix 
#' * y: The vector of outcomes
#' * std_X: standardized design matrix 
#' * std_X_details: a list with 2 vectors: 
#'    * 'center' (values used to center X)
#'    * 'scale' (values used to scale X)
#' * s: vector with the eigenvalues of K 
#' * U: the eigenvectors of K (same as left singular values of X). 
#' * ns: the indices for the nonsingular values of X
#' * penalty.factor: the penalty factors for the penalized non-singular values 
#' * snp_names: formatted column names of the design matrix 
#'
#' @keywords internal
#'
#' @examples
#' 
#' \dontrun{
#' # this is an internal function; to call this, you would need to use the triple 
#' # colon, e.g., plmmr:::plmm_prep()
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
#' prep2 <- plmm_prep(X = admix$X, y = admix$y, diag_K = TRUE, trace = TRUE)
#' }
#' 
plmm_prep <- function(std_X,
                      std_X_n,
                      std_X_p,
                      n,
                      p,
                      y,
                      k = NULL,
                      K = NULL,
                      diag_K = NULL,
                      eta_star = NULL,
                      fbm_flag,
                      penalty.factor = rep(1, ncol(X)),
                      trace = NULL, 
                      ...){
  
  
  ## coersion
  U <- s <- eta <- NULL
  
  if(is.null(k)){
    trunc_flag <- FALSE
  }
  
  if(!is.null(k)){
    
    if(k < min(std_X_n, std_X_p)){
      trunc_flag <- TRUE
    } else {
      stop("\nk value is out of bounds. 
         \nIf specified, k must be an integer in the range from 1 to min(nrow(X), ncol(X)). 
         \nwhere X does not include singular columns. For help detecting singularity,
         \nsee ncvreg::std()")
    }
    
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
    # TODO: add a check for list names
    # stopifnot("K_approx" %in% names(K))
    # sign_check <- flip_signs(X = K$K_approx, U = U, V = U, d = s)
    # U <- sign_check$U
  }
 
  # otherwise, need to do SVD:
  if(sum(c(flag1, flag2, flag3)) == 0){
    if(trace){cat("\nStarting decomposition.")}
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    if(is.null(K) & is.null(s)){
      # NB: the is.null(s) keeps you from overwriting the 3 preceding special cases 
      
      # approach to decomposition:
      # n > p: take SVD of X
      # n <= p: construct K, then take eigen(K)
      if(std_X_n > std_X_p){
      if(trace){cat("\nSince n > p, PLMM is calculating the SVD of X.")}
        svd_res <- svd_X(std_X = std_X, k = k, trunc_flag = trunc_flag,
                         fbm_flag = fbm_flag, trace = trace)
        s <- (svd_res$d^2)*(1/std_X_p)
        U <- svd_res$U
      
      } else if (std_X_n <= std_X_p){
        if(trace){cat("\nSince p > n, PLMM is calculating the eigendecomposition of K")}
        eigen_res <- eigen_K(std_X, p, fbm_flag = fbm_flag) 
        s <- eigen_res$s
        U <- eigen_res$U
        # check signs 
        # sign_check <- flip_signs(X = eigen_res$K, U = U, V = U, d = s)
        # U <- sign_check$U
      }
      
    } else {
      # last case: K is a user-supplied matrix 
      eigen_res <- eigen(K)
      s <- eigen_res$values*(1/std_X_p) # note: our definition of the RRM averages over the number of features
      U <- eigen_res$vectors
      # check signs
      # sign_check <- flip_signs(X = K, U = U, V = U, d = s)
      # U <- sign_check$U
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
    eta <- estimate_eta(n = std_X_n, s = s, U = U, y = y) 
  } else {
    # otherwise, use the user-supplied value (this option is mainly for simulation)
    eta <- eta_star
  }

# if FBM, keep U filebacked
  if(fbm_flag){
    U <- bigstatsr::as_FBM(U)
  }

  # return values to be passed into plmm_fit(): 
  ret <- structure(list(
    n = n,
    p = p,
    std_X_p = std_X_p,
    y = y,
    std_X = std_X,
    y = y,
    s = s,
    U = U,
    eta = eta, # carry eta over to fit 
    trace = trace))
  
 
  return(ret)
  
  
  
}
