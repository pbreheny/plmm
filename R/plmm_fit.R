#' PLMM fit: a function that fits a PLMM using the values returned by plmm_prep()
#' This is an internal function for \code{cv.plmm}
#' @param prep A list as returned from \code{plmm_prep}
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100. 
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @return A list with these components: 
#' * std_X: The standardized design matrix 
#' * SUX: first partial result of data rotation 
#' * SUy: second partial result of data rotation 
#' * eta: numeric value representing the ratio of variances. 
#' * std_SUX: re-standardized rotated design matrix. This is 'fed' into \code{plmm_fit()}. 
#' * b: The values returned in the 'beta' argument of the ncvfit() object
#' * lambda: The sequence of lambda values used in model fitting 
#' * iter: The number of iterations at each given lambda value 
#' * converged: The convergence status at each given lambda value 
#' * penalty: The type of penalty used in model fitting
#' * penalty.factor: A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' * ns: The indices of the non-singular columns of the ORIGINAL design matrix
#' * ncol_X: The number of columns in the ORIGINAL design matrix 
#' 
#' @keywords internal 
#'

plmm_fit <- function(prep, 
                     penalty = "MCP",
                     gamma,
                     alpha = 1,
                     # lambda.min = ifelse(n>p, 0.001, 0.05),
                     lambda.min,
                     nlambda = 100,
                     lambda,
                     eps = 1e-04,
                     max.iter = 10000,
                     convex = TRUE,
                     warn = TRUE,
                     init = NULL,
                     returnX = TRUE){
  
  # set default gamma (will need this for cv.plmm)
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  # set default init
  if(is.null(init)) init <- rep(0, prep$ncol_X)
  
  # error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(init)!=prep$ncol_X) stop("Dimensions of init and X do not match", call.=FALSE)
  
  if(prep$trace){cat("Beginning standardization + rotation.")}
  
  # estimate eta if needed
  if (is.null(prep$eta)) {
    eta <- estimate_eta(S = prep$S, U = prep$U, y = prep$y) 
  } else {
    # otherwise, use the user-supplied value (this is mainly for simulation)
    eta <- prep$eta
  }
  
  # rotate data
  W <- diag((eta * prep$S + (1 - eta))^(-1/2), nrow = length(prep$S)) 
  SUX <- W %*% crossprod(prep$U, cbind(1, prep$std_X)) # add column of 1s for intercept
  SUy <- drop(W %*% crossprod(prep$U, prep$y))
  
  # re-standardize rotated SUX
  std_SUX_temp <- scale_varp(SUX[,-1, drop = FALSE])
  std_SUX_noInt <- std_SUX_temp$scaled_X
  std_SUX <- cbind(SUX[,1, drop = FALSE], std_SUX_noInt) # re-attach intercept
  attr(std_SUX,'scale') <- std_SUX_temp$scale_vals
  
  # calculate population var without mean 0; will need this for call to ncvfit()
  xtx <- apply(std_SUX, 2, function(x) mean(x^2, na.rm = TRUE)) 
  
  if(prep$trace){cat("Setup complete. Beginning model fitting.\n")}
  
  # remove initial values for coefficients representing columns with singular values
  init <- init[prep$ns] 
  
  # set up lambda
  if (missing(lambda)) {
    lambda <- setup_lambda(X = std_SUX,
                           y = SUy,
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min = lambda.min,
                           penalty.factor = prep$penalty.factor)
    user.lambda <- FALSE
  } else {
    # make sure (if user-supplied sequence) is in DESCENDING order
    if(length(lambda) > 1){
      if (max(diff(lambda)) > 0) stop("User-supplied lambda sequence must be in descending (largest -> smallest) order")
    }
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  # make sure to *not* penalize the intercept term 
  new.penalty.factor <- c(0, prep$penalty.factor)
  
  # placeholders for results
  init <- c(0, init) # add initial value for intercept
  resid <- drop(SUy - std_SUX %*% init)
  b <- matrix(NA, nrow=ncol(std_SUX), ncol=nlambda) 
  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  
  # main attraction 
  ## set up progress bar -- this can take a while
  if(prep$trace){pb <- txtProgressBar(min = 0, max = nlambda, style = 3)}
  ## TODO: think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(std_SUX, SUy, init, resid, xtx, penalty, gamma, alpha, lam, eps, max.iter, new.penalty.factor, warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    resid <- res$resid
    if(prep$trace){setTxtProgressBar(pb, ll)}
  }
  
  # reconstruct K to calculate V 
  # this is on the standardized X scale 
  estimated_V <- eta * tcrossprod(prep$U %*% diag(prep$S), prep$U) + (1-eta)*diag(nrow = nrow(prep$U)) 
  
  ret <- structure(list(
    std_X = prep$std_X,
    y = prep$y,
    S = prep$S,
    U = prep$U,
    std_SUX = std_SUX,
    SUX = SUX,
    SUy = SUy,
    S = prep$S,
    U = prep$U,
    ncol_X = prep$ncol_X, 
    nrow_X = prep$nrow_X, 
    lambda = lambda,
    b = b,
    eta = eta,
    iter = iter,
    converged = converged, 
    loss = loss, 
    penalty = penalty, 
    penalty.factor = new.penalty.factor,
    gamma = gamma,
    alpha = alpha,
    ns = prep$ns,
    snp_names = prep$snp_names,
    penalty = penalty,
    nlambda = nlambda,
    eps = eps,
    max.iter = max.iter,
    warn = warn,
    init = init,
    returnX = prep$returnX,
    trace = prep$trace, 
    estimated_V = estimated_V
  )) 
  
  return(ret)
  
  
  
}
