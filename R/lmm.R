#' a function to implement a mixed model for structured data *without* penalization
#' 
#' NB: this function is simply a wrapper for lmm_prep -> lmm_fit -> lmm_format
#' @param X Design matrix. MUST NOT be high dimensional -- the modeling here will call `lm()`. For high-dimensional data, use `plmm()` instead.
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list as created by choose_k.
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @return A list including the estimated coefficients on the original scale, as well as other model fitting details 
#' 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' # compare basic LM with LMM (latter accounts for family relationships)
#' fit1 <- lmm(X = pedigree$X, y = pedigree$clinical$y, diag_K = TRUE)
#' fit2 <- lmm(X = pedigree$X, y = pedigree$clinical$y, K = pedigree$K)
#' fit1$beta_vals; fit2$beta_vals

lmm <- function(X,
                 y,
                 k = NULL, 
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 eps = 1e-04,
                 max.iter = 10000,
                 dfmax = ncol(X) + 1,
                 warn = TRUE,
                 init = rep(0, ncol(X)),
                 returnX = TRUE,
                 trace = FALSE) {
  
  ## check types 
  if ("SnpMatrix" %in% class(X)) X <- methods::as(X, 'numeric')
  if("FBM.code256" %in% class(X)) stop("plmm does not work with FBM objects at this time. This option is in progress. For now, X must be a numeric matrix.")
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
  
  # working with user-specified K
  if (!is.null(K)){
    # first, check type/class:
    if (!inherits(K, "matrix") & !is.list(K)) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be either (1) able to be coerced to a matrix or (2) be a list.", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(X) <- "double" # change K to X 
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (is.list(K)) {
      if(!('d' %in% names(K) & 'U' %in% names(K))){stop('Components d and U not both found in list supplied for K.')}
    }
    # last thing: check dimensions
    if (is.matrix(K)){
      if (dim(K)[1] != nrow(X) || dim(K)[2] != nrow(X)) stop("Dimensions of K and X do not match", call.=FALSE)
    } else if (is.list(K)) {
      if (nrow(X) != nrow(K$U)) stop("Dimensions of K and X do not match.")
    }
    
  }
  
  if(trace){cat("Passed all checks. Beginning singular value decomposition.\n")}
  
  the_prep <- lmm_prep(X = X,
                        y = y,
                        K = K,
                        k = k,
                        diag_K = diag_K,
                        eta_star = eta_star,
                        returnX = returnX,
                        trace = trace)
  
  
  if(trace){cat("Beginning model fitting.\n")}
  the_fit <- lmm_fit(prep = the_prep,
                      eps = eps,
                      max.iter = max.iter,
                      warn = warn,
                      init = init,
                      returnX = returnX)
  
  if(trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  the_final_product <- lmm_format(fit = the_fit, 
                                   dfmax = dfmax, 
                                   X = X)
  
  return(the_final_product)
  
  
}
