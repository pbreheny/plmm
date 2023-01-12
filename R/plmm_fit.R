#' PLMM fit: a function that fits a PLMM using the values returned by plmm_prep()
#' This is an internal function for \code{cv.plmm}
#' @param X The same design matrix used in \code{plmm_prep}
#' @param y The continuous outcome vector used in \code{plmm_prep}
#' @param prep A list as returned from \code{plmm_prep}
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
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @return
#' @export
#'
#' @examples
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' fit1 <- plmm_fit(X = admix$X, y = admix$y, prep = prep1)
plmm_fit <- function(X,
                     y,
                     prep, 
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
                     init = rep(0, ncol(X)),
                     returnX = TRUE, 
                     trace = FALSE){
  
  # coercion
  penalty <- match.arg(penalty)
  
  # set default gamma
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  # error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(init)!=ncol(X)) stop("Dimensions of init and X do not match", call.=FALSE)
  
  # designate the dimensions of the design matrix 
  p <- ncol(X) 
  n <- nrow(X)
  
  # re-standardize rotated SUX
  std_SUX_temp <- scale_varp(prep$SUX[,-1, drop = FALSE])
  std_SUX_noInt <- std_SUX_temp$scaled_X
  std_SUX <- cbind(prep$SUX[,1, drop = FALSE], std_SUX_noInt) # re-attach intercept
  attr(std_SUX,'scale') <- std_SUX_temp$scale_vals
  
  # calculate population var without mean 0; will need this for call to ncvfit()
  xtx <- apply(std_SUX, 2, function(x) mean(x^2, na.rm = TRUE)) 
  
  if(trace){cat("Setup complete. Beginning model fitting.\n")}
  
  # remove initial values for coefficients representing columns with singular values
  init <- init[prep$ns] 
  
  # set up lambda
  if (missing(lambda)) {
    lambda <- setup_lambda(X = std_SUX,
                           y = prep$SUy,
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
  resid <- drop(prep$SUy - std_SUX %*% init)
  b <- matrix(NA, nrow=ncol(std_SUX), ncol=nlambda) 
  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  
  # main attraction 
  ## set up progress bar -- this can take a while
  if(trace){pb <- txtProgressBar(min = 0, max = nlambda, style = 3)}
  ## TODO: think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(std_SUX, prep$SUy, init, resid, xtx, penalty, gamma, alpha, lam, eps, max.iter, new.penalty.factor, warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    resid <- res$resid
    if(trace){setTxtProgressBar(pb, ll)}
  }
  
  # eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")
  convex.min <- if (convex) convexMin(b, std_SUX, penalty, gamma, lambda*(1-alpha), family = 'gaussian', new.penalty.factor) else NULL
  
  # reverse the transformations of the beta values 
  beta_vals <- untransform(b, prep$ns, X, prep$std_X, prep$SUX, std_SUX)
  
  if(trace){cat("\nBeta values are estimated -- almost done!\n")}
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows, lambda values on the columns
  varnames <- if (is.null(colnames(X))) paste("K", 1:ncol(X), sep="") else colnames(X)
  varnames <- c("(Intercept)", varnames)
  dimnames(beta_vals) <- list(varnames, lamNames(lambda))
  
  ## output
  val <- structure(list(beta_vals = beta_vals,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        new.penalty.factor = new.penalty.factor,
                        ns_idx = c(1, 1 + prep$ns), # PAY ATTENTION HERE! 
                        iter = iter,
                        converged = converged),
                   class = "plmm")
  if (missing(returnX)) {
    if (utils::object.size(prep$SUX) > 1e8) {
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
  } 
  return(val)
  
  
}