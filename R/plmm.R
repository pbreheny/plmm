#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' NB: this function is simply a wrapper for plmm_prep -> plmm_fit -> plmm_format
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p)`
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}).
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param k An integer between 1 and \code{nrow(K)} indicating the number of singular values requested *if* package \code{RSpectra} is installed. Defaults to NULL. 
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/Spenncath penalty (see details). Default is 3 for MCP and 3.7 for Spenncath.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/Spenncath penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/Spenncath penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100. 
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param svd_details Logical: should the details from the SVD, such as the singular values (S) and singular vectors (U) be returned? Default is TRUE. 
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
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
#' # using admix data 
#' fit_admix1 <- plmm(X = admix$X, y = admix$y)
#' s1 <- summary(fit_admix1, idx = 99)
#' print(s1)
#' 
#' # using admix data and k = 50 
#' fit_admix2 <- plmm(X = admix$X, y = admix$y, k = 50)
#' s2 <- summary(fit_admix2, idx = 99)
#' print(s2)
#' 
#' # now use PLINK data files
#' \dontrun{
#' 
#' penncath_mid <- process_plink(prefix = "penncath_mid", dataDir = plink_example(path="penncath_mid.fam", parent=T))
#' penncath_clinical <- read.csv(plink_example(path="penncath_clinical.csv"))
#' # for the sake of illustration, I use a simple mean imputation for the outcome 
#' penncath_clinical$hdl_impute <- ifelse(is.na(penncath_clinical$hdl), mean(penncath_clinical$hdl, na.rm = T), penncath_clinical$hdl)
#' 
#' # fit with no 'k' specified
#' fit_plink1 <- plmm(X = penncath_mid$genotypes, y = penncath_clinical$hdl_impute, trace = TRUE)
#' summary(fit_plink1, idx = 5)
#' # Runs in ~219 seconds (3.65 mins) on my 2015 MacBook Pro
#' 
#' # fit with 'k = 5' specified (so using RSpectra::svds())
#' fit_plink2 <- plmm(X = penncath_mid$genotypes, y = penncath_clinical$hdl_impute, k = 5, trace = TRUE)
#' # Runs in ~44 seconds on my 2015 MacBook Pro
#' summary(fit_plink2, idx = 5);summary(fit_plink2, idx = 95)
#' }


plmm <- function(X,
                 y,
                 K = NULL,
                 k = min(nrow(X), ncol(X)), 
                 # TODO: think more deeply about what the best choice is for the 
                 # default value of k
                 eta_star = NULL,
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma,
                 alpha = 1,
                 # lambda.min = ifelse(n>p, 0.001, 0.05),
                 lambda.min,
                 nlambda = 100,
                 lambda,
                 svd_details = TRUE,
                 eps = 1e-04,
                 max.iter = 10000,
                 convex = TRUE,
                 dfmax = ncol(X) + 1,
                 warn = TRUE,
                 penalty.factor = rep(1, ncol(X)),
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
  
  # coercion for penalty
  penalty <- match.arg(penalty)
  
  # set default gamma
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  
  if(trace){cat("Passed all checks. Beginning singular value decomposition.\n")}
  
  the_prep <- plmm_prep(X = X,
                        y = y,
                        K = K,
                        k = k,
                        eta_star = eta_star,
                        penalty.factor = penalty.factor,
                        returnX = returnX,
                        trace = trace)
  
  
  if(trace){cat("Beginning model fitting.\n")}
  the_fit <- plmm_fit(prep = the_prep,
                      penalty = penalty,
                      gamma = gamma,
                      alpha = alpha,
                      lambda.min = lambda.min,
                      nlambda = nlambda,
                      lambda = lambda,
                      svd_details = svd_details,
                      eps = eps,
                      max.iter = max.iter,
                      warn = warn,
                      init = init,
                      returnX = returnX)
  
  if(trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  the_final_product <- plmm_format(fit = the_fit, 
                                   convex = convex,
                                   dfmax = dfmax, 
                                   X = X, K = K)
  
  return(the_final_product)
  
  
}
