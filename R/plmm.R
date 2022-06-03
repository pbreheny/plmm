#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param V Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (typically computed as XX^T).
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If V is a known covariance matrix that is full rank, this should be 1.
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
#' @param init Initial values for coefficients. Default is TRUE. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @param intercept Logical flag for whether an intercept should be included.
#' @param standardizeX Logical flag for X matrix standardization, prior to data transformation. The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize.
#' @param standardizeRtX Logical flag for transformed X matrix scaling (or rescaling). The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize, but this is generally recommended.
#' @param rotation Logical flag to indicate whether the weighted rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @param ... Not used.
#' @importFrom zeallot %<-%
#' @export


plmm <- function(X,
                 y,
                 V,
                 eta_star,
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma,
                 alpha = 1,
                 lambda.min = ifelse(n>p, 0.001, 0.05),
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max.iter = 10000,
                 convex = TRUE,
                 dfmax = p + 1,
                 warn = TRUE,
                 penalty.factor = rep(1, ncol(X)),
                 init = rep(0, ncol(X)),
                 returnX = TRUE,
                 intercept = TRUE,
                 standardizeX = TRUE,
                 standardizeRtX = TRUE,
                 rotation = TRUE,
                 ...) {

  # Coersion
  U <- S <- SUX <- SUy <- eta <- NULL
  penalty <- match.arg(penalty)
  if (missing(V)) stop('Similarity matrix must be provided.')
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
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

  # Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  if (length(init)!=ncol(X)) stop("Dimensions of init and X do not match", call.=FALSE)
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  if (!missing(V)){
    if (!inherits(V, "matrix")) {
      tmp <- try(V <- stats::model.matrix(~0+., data=V), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("V must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(V)=="integer") storage.mode(X) <- "double"
    if (typeof(V)=="character") stop("V must be a numeric matrix", call.=FALSE)
    if (dim(V)[1] != nrow(X) || dim(V)[2] != nrow(X)) stop("Dimensions of V and X do not match", call.=FALSE)
  }

  ## Deprecation support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Standardize X
  if (standardizeX){
    XX <- ncvreg::std(X)
  } else {
    XX <- X
    attributes(XX)$center <- rep(0, ncol(XX))
    attributes(XX)$scale <- rep(1, ncol(XX))
    attributes(XX)$nonsingular <- 1:ncol(XX)
    # X_ns <- attr(XX, "nonsingular") # shouldn't need this value in this scenario
  }
  ns <- attr(XX, "nonsingular")
  init <- init[ns] # remove any singular values

  penalty.factor <- penalty.factor[ns]

  p <- ncol(XX)
  n <- nrow(XX)

  ## Rotate data
  if (!missing(eta_star)){
    c(SUX, SUy, eta, U, S) %<-% rotate_data(XX, y, V, intercept, rotation, eta_star)
  } else {
    c(SUX, SUy, eta, U, S) %<-% rotate_data(XX, y, V, intercept, rotation)
  }

  if (intercept) penalty.factor <- c(0, penalty.factor)

  ## Re-standardize rotated SUX
  if (standardizeRtX){
    if (intercept){
      SUXX_noInt <- scale_varp(SUX[,-1, drop = FALSE])
      attributes(SUXX_noInt)$nonsingular <- attributes(ncvreg::std(SUX[,-1, drop = FALSE]))$nonsingular
      SUXX <- cbind(SUX[,1, drop = FALSE], SUXX_noInt)
    } else {
      SUXX_noInt <- scale_varp(SUX)
      attributes(SUXX_noInt)$nonsingular <- attributes(ncvreg::std(SUX))$nonsingular
      SUXX <- SUXX_noInt
    }
    attributes(SUXX)$scale <- attr(SUXX_noInt, 'scale')
    attributes(SUXX)$nonsingular <- attr(SUXX_noInt, 'nonsingular')
  } else {
    SUXX <- SUX
    if (intercept){
      attributes(SUXX)$scale <- rep(1, ncol(SUX) - 1)
      attributes(SUXX)$nonsingular <- 1:(ncol(SUX) - 1)
      # xtx <- c(1, apply(SUXX[,-1], 2, varp))
    } else {
      attributes(SUXX)$scale <- rep(1, ncol(SUX))
      attributes(SUXX)$nonsingular <- 1:ncol(SUX)
    }
  }
  xtx <- apply(SUXX, 2, function(x) mean(x^2, na.rm = TRUE)) # population var without mean 0

  ## Set up lambda
  if (missing(lambda)) {
    lambda <- setup_lambda(SUXX, SUy, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Placeholders for results

  if (intercept) init <- c(0, init) # add initial value for intercept
  resid <- drop(SUy - SUXX %*% init)
  b <- matrix(NA, ncol(SUXX), nlambda)
  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  # think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(SUXX, SUy, init, resid, xtx, penalty, gamma, alpha, lam, eps, max.iter, penalty.factor, warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    resid <- res$resid
  }

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")
  convex.min <- if (convex) convexMin(b, SUXX, penalty, gamma, lambda*(1-alpha), family = 'gaussian', penalty.factor) else NULL

  # unscale transformed data
  bb <- unscale(b[, ind, drop = FALSE], SUX, SUXX, intercept)

  # unstandardize original data
  beta <- unstandardize(bb, X, XX, intercept)

  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  if (intercept) varnames <- c("(Intercept)", varnames)
  dimnames(beta) <- list(varnames, lamNames(lambda))

  ## Output
  val <- structure(list(beta = beta,
                        eta = eta,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n,
                        iter = iter,
                        converged = converged),
                        class = "plmm")
  if (missing(returnX)) {
    if (utils::object.size(SUX) > 1e8) {
      warning("Due to the large size of SUX (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$SUX <- SUX
    val$SUy <- SUy
    val$U <- U
    val$S <- S
    val$XX <- XX
  }
  return(val)
}
