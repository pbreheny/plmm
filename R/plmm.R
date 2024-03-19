#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' NB: this function is simply a wrapper for plmm_prep -> plmm_fit -> plmm_format
#' @param X Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to `get_data()`. 
#' * Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
#' @param y Continuous outcome vector. Logistic regression modeling is still in development.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for Spenncath.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/Spenncath penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/Spenncath penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100. 
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @return A list which includes: 
#'  * `beta_vals`: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * `rotated_scale_beta_vals`: the matrix of estimated coefficients on the ~rotated~ scale. This is the scale on which the model was fit. 
#'  * `lambda`: a numeric vector of the lasso tuning parameter values used in model fitting. 
#'  * `eta`: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
#'  * `s`: a vectof of the eigenvalues of relatedness matrix `K`; see `relatedness_mat()` for details.
#'  * `U`: a matrix of the eigenvalues of relatedness matrix `K`
#'  * `rot_y`: the vector of outcome values on the rotated scale. This is the scale on which the model was fit. 
#'  * `linear.predictors`: the matrix resulting from the product of `stdrot_X` and the estimated coefficients on the ~rotated~ scale.
#'  * `penalty`: character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * `gamma`: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
#'  * `alpha`: numeric value indicating the elastic net tuning parameter. 
#'  * `convex.min`: NULL (This is an option we will add in the future!)
#'  * `loss`: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `penalty.factor`: vector of indicators corresponding to each predictor, where 1 = predictor was penalized. 
#'  * `ns_idx`: vector with the indicies of predictors which were constant features (i.e., had no variation).
#'  * `p`: the number of features 
#'  * `n`: the number of observations (instances)
#'  * `iter`: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: vector of logical values indicating whether the model fitting converged at each value of `lambda`
#' 
#' @export
#' 
#' @examples 
#' # using admix data 
#' fit_admix1 <- plmm(X = admix$X, y = admix$y)
#' s1 <- summary(fit_admix1, idx = 99)
#' print(s1)
#' plot(fit_admix1)
#' 
#' # an example with p > n:
#' fit_admix2 <- plmm(X = admix$X[1:50, ], y = admix$y[1:50])
#' s2 <- summary(fit_admix2, idx = 99)
#' print(s2)
#' plot(fit_admix2) # notice: the default penalty is MCP
#' 
plmm <- function(X,
                 y,
                 k = NULL, 
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma,
                 alpha = 1,
                 # lambda.min = ifelse(n>p, 0.001, 0.05),
                 lambda.min,
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max.iter = 10000,
                 convex = TRUE,
                 dfmax = ncol(X) + 1,
                 warn = TRUE,
                 penalty.factor = rep(1, ncol(X)),
                 init = rep(0, ncol(X)),
                 trace = FALSE) {

  ## check types 
  if("character" %in% class(X)){
    dat <- get_data(path = X, returnX = TRUE, trace = trace)
    X <- dat$X
  }
  if ("SnpMatrix" %in% class(X)) X <- methods::as(X, 'numeric')
  if("FBM.code256" %in% class(X)) stop("plmm does not work with FBM objects at this time. This option is in progress. \nFor now, design matrix X must be a numeric matrix.")
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
  
  if (!is.null(K)){
    # first, check type/class:
    if (!inherits(K, "matrix") & !is.list(K)) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be either (1) able to be coerced to a matrix or (2) be a list with elements 's' and 'U'.", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(X) <- "double" # change K to X 
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (is.list(K)) {
      if(!('s' %in% names(K) & 'U' %in% names(K))){stop('Components s and U not both found in list supplied for K.')}
    }
    # last thing: check dimensions
    if (is.matrix(K)){
      if (dim(K)[1] != nrow(X) || dim(K)[2] != nrow(X)) stop("Dimensions of K and X do not match", call.=FALSE)
    } else if (is.list(K)) {
      if (nrow(X) != nrow(K$U)) stop("Dimensions of K and X do not match.")
    }
  
  }
  # warn about computational time for large K 
  if(is.null(k) & is.null(diag_K) & (nrow(X) > 1000)){
    warning("The number of observations is large, and k is not specified.
    \nThis can dramatically increase computational time -- the SVD calculation is expensive.
            \nIf the observations are unrelated, please set diag_K = TRUE. SVD is not needed in this case.
            \nOtherwise, consider using choose_k() first to get an approximation for your relatedness matrix.")
  }
  # coercion for penalty
  penalty <- match.arg(penalty)
  
  # set default gamma
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  
  if(trace){cat("\nPassed all checks. Beginning singular value decomposition.\n")}
  the_prep <- plmm_prep(X = X,
                        y = y,
                        K = K,
                        k = k,
                        diag_K = diag_K,
                        eta_star = eta_star,
                        penalty.factor = penalty.factor,
                        trace = trace)

  if(trace){cat("Beginning model fitting.\n")}
  the_fit <- plmm_fit(prep = the_prep,
                      penalty = penalty,
                      gamma = gamma,
                      alpha = alpha,
                      lambda.min = lambda.min,
                      nlambda = nlambda,
                      lambda = lambda,
                      eps = eps,
                      max.iter = max.iter,
                      warn = warn,
                      init = init,
                      convex = convex,
                      dfmax = dfmax)

  if(trace){cat("\nBeta values are estimated -- almost done!")}
  the_final_product <- plmm_format(fit = the_fit, X = X)
  
  return(the_final_product)
  
  
}
