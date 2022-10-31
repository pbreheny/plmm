#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}).
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param k An integer between 1 and \code{nrow(K)} indicating the number of singular values requested *if* package \code{RSpectra} is installed. Defaults to NULL. 
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
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
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' 
#' @return A list including the estimated coeficients on the original scale, as well as other model fitting details 
#' 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' s <- summary(fit, idx = 99)
#' print(s)
#' 
#' # now use PLINK data files
#' \dontrun{
#' cad <- process_plink(prefix = "cad", dataDir = plink_example(path="cad.fam", parent=T))
#' cad_clinical <- read.csv((plink_example(path="cad_clinical.csv"))
#' # for the sake of illustration, I use a simple mean imputation for the outcome 
#' cad_clinical$hdl_impute <- ifelse(is.na(cad_clinical$hdl), mean(cad_clinical$hdl, na.rm = T), cad_clincal$hdl)
#' fit_plink <- plmm(X = coerce_snpmatrix(cad$genotypes), y = cad_clinical$hdl_impute, k = 10)
#' }


plmm <- function(X,
                 y,
                 K,
                 k = NULL,
                 eta_star,
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
                 dfmax = p + 1,
                 warn = TRUE,
                 penalty.factor = rep(1, ncol(X)),
                 init = rep(0, ncol(X)),
                 returnX = TRUE) {

  ## coersion
  U <- S <- SUX <- SUy <- eta <- NULL
  penalty <- match.arg(penalty)
  
  ## set defaults 
  if(missing(K)){K <- relatedness_mat(X)}  
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
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
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)
  
  ## error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  if (length(init)!=ncol(X)) stop("Dimensions of init and X do not match", call.=FALSE)
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  if (!missing(K)){
    if (!inherits(K, "matrix")) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(X) <- "double" # change K to X 
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (dim(K)[1] != nrow(X) || dim(K)[2] != nrow(X)) stop("Dimensions of K and X do not match", call.=FALSE)
  }

## standardize X
  # NB: the following line will eliminate singular columns (eg monomorphic SNPs)
  #  from the design matrix. 
  std_X <- ncvreg::std(X)
  
  
  # identify nonsingular values in the standardized X matrix  
  ns <- attr(std_X, "nonsingular")
  
  # remove initial values for coefficients representing columns with singular values
  init <- init[ns] 

## keep only those penalty factors which penalize non-singular values 
  penalty.factor <- penalty.factor[ns]

## designate the dimensions of the design matrix 
  p <- ncol(std_X) 
  n <- nrow(std_X)

## rotate data
  if (!missing(eta_star)){
    if(is.null(k)){
      c(SUX, SUy, eta, U, S) %<-% rotate_data(std_X, y, K, eta_star)
    } else {
      c(SUX, SUy, eta, U, S) %<-% rotate_data(std_X, y, K, eta_star, k)
    }
    
  } else {
    if(is.null(k)){
      c(SUX, SUy, eta, U, S) %<-% rotate_data(std_X, y, K)
    } else {
      c(SUX, SUy, eta, U, S) %<-% rotate_data(std_X, y, K, k)
    }
    
  }

## re-standardize rotated SUX
std_SUX_temp <- scale_varp(SUX[,-1, drop = FALSE])
std_SUX_noInt <- std_SUX_temp$scaled_X

std_SUX <- cbind(SUX[,1, drop = FALSE], std_SUX_noInt) # re-attach intercept
attr(std_SUX,'scale') <- std_SUX_temp$scale_vals

## calculate population var without mean 0; will need this for call to ncvfit()
xtx <- apply(std_SUX, 2, function(x) mean(x^2, na.rm = TRUE)) 

## set up lambda
if (missing(lambda)) {
    lambda <- setup_lambda(X = std_SUX,
                           y = SUy,
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min = lambda.min,
                           penalty.factor = penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

# make sure lambda sequence (if user-supplied) is in DESCENDING order
if (max(diff(lambda)) > 0) stop("User-supplied lambda sequence must be in descending (largest -> smallest) order")

  
# make sure to *not* penalize the intercept term 
penalty.factor <- c(0, penalty.factor)
  

## placeholders for results
init <- c(0, init) # add initial value for intercept
resid <- drop(SUy - std_SUX %*% init)
b <- matrix(NA, nrow=ncol(std_SUX), ncol=nlambda) 
iter <- integer(nlambda)
converged <- logical(nlambda)
loss <- numeric(nlambda)

## main attraction 
  # think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(std_SUX, SUy, init, resid, xtx, penalty, gamma, alpha, lam, eps, max.iter, penalty.factor, warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    resid <- res$resid
  }

## eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")
  convex.min <- if (convex) convexMin(b, std_SUX, penalty, gamma, lambda*(1-alpha), family = 'gaussian', penalty.factor) else NULL

# reverse the transformations of the beta values 
beta_vals <- untransform(b, ns, X, std_X, SUX, std_SUX)

# give the matrix of beta_values readable names 
# SNPs (or covariates) on the rows, lambda values on the columns
varnames <- if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)
varnames <- c("(Intercept)", varnames)
dimnames(beta_vals) <- list(varnames, lamNames(lambda))

## output
val <- structure(list(beta_vals = beta_vals,
                        eta = eta,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n,
                        ns_idx = c(1, 1 + ns), # PAY ATTENTION HERE! 
                        iter = iter,
                        converged = converged),
                        class = "plmm")
  if (missing(returnX)) {
    if (utils::object.size(SUX) > 1e8) {
      warning("Due to the large size of SUX (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      # if it fits, it ships 
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- X # this is the original design matrix WITHOUT the intercept!
    val$y <- y
    
    val$std_X <- std_X
    
    val$U <- U
    val$S <- S
    
    val$SUX <- SUX
    val$SUy <- SUy
    
  } 
  return(val)
}
