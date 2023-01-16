#' PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model
#' This is an internal function for \code{cv.plmm}
#'
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}).
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
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
#' * penalty.factor: the penalty factors for the penalized non-singular values 
#' * snp_names: Formatted column names of the design matrix 
#' @export
#'
#' @examples
#' 
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
plmm_prep <- function(X,
                      y,
                      K = NULL,
                      eta_star,
                      penalty.factor = rep(1, ncol(X)),
                      returnX = TRUE,
                      trace = FALSE, ...){
  
  
  ## coersion
  U <- S <- SUX <- SUy <- eta <- NULL
  
  
  ## check types 
  if ("SnpMatrix" %in% class(X)) X <- methods::as(X, 'numeric')
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }

  # error checking 
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  
  # working with user-specified K
  if (!is.null(K)){
    if (!inherits(K, "matrix")) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(X) <- "double" # change K to X 
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (dim(K)[1] != nrow(X) || dim(K)[2] != nrow(X)) stop("Dimensions of K and X do not match", call.=FALSE)
  }
  
  if(trace){cat("Passed all checks. Beginning singular value decomposition.\n")}
  
  # standardize X
  # NB: the following line will eliminate singular columns (eg monomorphic SNPs)
  #  from the design matrix. 
  std_X <- ncvreg::std(X) 
  
  # identify nonsingular values in the standardized X matrix  
  ns <- attr(std_X, "nonsingular")
  
  # keep only those penalty factors which penalize non-singular values 
  penalty.factor <- penalty.factor[ns]
  
  # designate the dimensions of the design matrix 
  p <- ncol(X) 
  n <- nrow(X)
  
  # calculate SVD
  ## case 1: K is not specified (default to realized relatedness matrix)
  if (is.null(K)){
    if(trace){cat("No K specified - will use default definition of the \n realized relatedness matrix.\n")}
    c(D, U) %<-% svd(std_X, nv = 0) # don't need V
    # TODO: add RSpectra option here (or, use bigsnpr SVD method)
    S <- (D^2)/p # singular values of K, the realized relationship matrix
    
  } else {
    ## case 2: K is user-specified 
    S <- U <- NULL
    c(S, U) %<-% svd(K, nv = 0) # again, don't need V 
  }

  
  # return values to be passed into plmm_fit(): 
  ret <- structure(list(ncol_X = ncol(X),
                        y = y,
                        std_X = std_X,
                        S = S,
                        U = U,
                        ns = ns,
                        penalty.factor = penalty.factor,
                        trace = trace,
                        returnX = returnX,
                        snp_names = if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)))
  
  return(ret)
  
  
  
}
