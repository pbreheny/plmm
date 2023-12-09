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
#' * rot_X: first partial result of data rotation 
#' * rot_y: second partial result of data rotation 
#' * eta: numeric value representing the ratio of variances. 
#' * std_rot_X: re-standardized rotated design matrix. This is 'fed' into \code{plmm_fit()}. 
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
  if(is.null(init)) init <- rep(0, prep$p)
  
  # error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(init)!=prep$p) stop("Dimensions of init and X do not match", call.=FALSE)
  
  if(prep$trace){cat("\nBeginning standardization + rotation.")}
  
  # estimate eta if needed
  if (is.null(prep$eta)) {
    eta <- estimate_eta(s = prep$s, U = prep$U, y = prep$y) 
  } else {
    # otherwise, use the user-supplied value (this is mainly for simulation)
    eta <- prep$eta
  }

  # rotate data
  w <- (eta * prep$s + (1 - eta))^(-1/2)
  wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")
  rot_X <- wUt %*% cbind(1, prep$std_X)
  rot_y <- wUt %*% prep$y
  # re-standardize rotated rot_X
  stdrot_X_temp <- scale_varp(rot_X[,-1, drop = FALSE])
  stdrot_X_noInt <- stdrot_X_temp$scaled_X
  stdrot_X <- cbind(rot_X[,1, drop = FALSE], stdrot_X_noInt) # re-attach intercept
  attr(stdrot_X,'scale') <- stdrot_X_temp$scale_vals

  # calculate population var without mean 0; will need this for call to ncvfit()
  xtx <- apply(stdrot_X, 2, function(x) mean(x^2, na.rm = TRUE)) 
  
  if(prep$trace){cat("\nSetup complete. Beginning model fitting.")}
  
  # remove initial values for coefficients representing columns with singular values
  init <- init[prep$ns] 

  # set up lambda
  if (missing(lambda)) {
    lambda <- setup_lambda(X = stdrot_X,
                           y = rot_y,
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min = lambda.min,
                           penalty.factor = prep$penalty.factor)
    user.lambda <- FALSE
  } else {
    # make sure (if user-supplied sequence) is in DESCENDING order
    if(length(lambda) > 1){
      if (max(diff(lambda)) > 0) stop("\nUser-supplied lambda sequence must be in descending (largest -> smallest) order")
    }
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  # make sure to *not* penalize the intercept term 
  new.penalty.factor <- c(0, prep$penalty.factor)
  
  # placeholders for results
  init <- c(0, init) # add initial value for intercept
  resid <- drop(rot_y - stdrot_X %*% init)
  linear.predictors <- matrix(NA, nrow = nrow(stdrot_X), ncol=nlambda)
  b <- matrix(NA, nrow=ncol(stdrot_X), ncol=nlambda) 
  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  
  # main attraction 
  ## set up progress bar -- this can take a while
  if(prep$trace){pb <- txtProgressBar(min = 0, max = nlambda, style = 3)}
  ## TODO: think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(stdrot_X, rot_y, init, resid, xtx, penalty, gamma, alpha, lam, eps, max.iter, new.penalty.factor, warn)
    b[, ll] <- init <- res$beta
    linear.predictors[,ll] <- stdrot_X%*%(res$beta)
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    resid <- res$resid
    if(prep$trace){setTxtProgressBar(pb, ll)}
  }
  ret <- structure(list(
    y = prep$y,
    p = prep$p, 
    n = prep$n, 
    s = prep$s,
    U = prep$U,
    rot_X = rot_X,
    rot_y = rot_y,
    stdrot_X = stdrot_X,
    lambda = lambda,
    b = b,
    linear.predictors = linear.predictors,
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
    trace = prep$trace)) 
  
  return(ret)
  
  
  
}
