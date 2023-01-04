#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized penalized linear mixed models over a grid of values for the regularization parameter lambda.
#' @param X Design matrix for model fitting. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector for model fitting.
#' @param K Known or estimated similarity matrix.
#' @param eta_star Optional arg. to \code{plmm_prep}. Defaults to NULL.
#' @param penalty.factor Optional arg. to \code{plmm_prep}. Defaults to 1 for all predictors (except the intercept). 
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'response'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'blup'} predictions are based on the linear predictor plus the estimated random effect (BLUP). Defaults to 'response'.
#' @param ... Additional arguments to \code{plmm_fit}
#' @param cluster cv.plmm can be run in parallel across a cluster using the parallel package. The cluster must be set up in advance using the makeCluster function from that package. The cluster must then be passed to cv.plmm.
#' @param nfolds The number of cross-validation folds. Default is 10.
#' @param fold Which fold each observation belongs to. By default the observations are randomly assigned.
#' @param seed You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY Should cv.plmm return the linear predictors from the cross-validation folds? Default is FALSE; if TRUE, this will return a matrix in which the element for row i, column j is the fitted value for observation i from the fold in which observation i was excluded from the fit, at the jth value of lambda.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @export
#' 
#' @examples 
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' cv_s <- summary.cv.plmm(cv_fit, lambda = "1se")
#' print(cv_s)


cv.plmm <- function(X,
                    y,
                    K,
                    eta_star = NULL,
                    penalty.factor = rep(1, ncol(X)),
                    type = 'response',
                    ...,
                    cluster,
                    nfolds=10,
                    seed,
                    fold,
                    returnY=FALSE,
                    trace=FALSE) {

  # Default type is 'response'
  if(missing(type)) {type == 'response'} 

  prep.args <- c(list(X = X,
                     y = y,
                     K = K,
                     eta_star = eta_star,
                     penalty.factor = penalty.factor))
  
  # implement preparation steps for model fitting 
  prep <- do.call('plmm_prep', prep.args)
  
  
  # implement full model fit 
  fit.args <- c(list(X = X, y = y, prep = prep), list(...))
  fit <- do.call('plmm_fit', fit.args)
  
  # set up arguments for cv 
  cv.args <- fit.args
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$lambda <- fit$lambda
  cv.args$eta_star <- fit$eta
  
  if (type == 'blup') {cv.args$returnX <- TRUE}
  
  n <- length(fit$y)
  E <- Y <- matrix(NA, nrow=nrow(X), ncol=length(fit$lambda))

  if (!missing(seed)) set.seed(seed)
  sde <- sqrt(.Machine$double.eps)
  if (missing(fold)) {
    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "K", "fold", "type", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(penalizedLMM))
    fold.results <- parallel::parLapply(cl=cluster, X=1:max(fold), fun=cvf, XX=X, y=y, fold=fold, type=type, cv.args=cv.args)
  }

  if (trace) cat("\nStarting cross validation\n")  
  # set up progress bar -- this can take a while
  if(trace){pb <- txtProgressBar(min = 0, max = nfolds, style = 3)}
  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
      if (trace) {setTxtProgressBar(pb, i)}
    } else {
      res <- cvf(i = i, XX = X, y = y, fold = fold, type = type, cv.args = cv.args)
      if (trace) {setTxtProgressBar(pb, i)}
    }
    # browser()
    E[fold==i, 1:res$nl] <- res$loss
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  ## Return min lambda idx
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(n)
  min <- which.min(cve)

  ## Return lambda 1se idx
  l.se <- cve[min] - cvse[min]
  u.se <- cve[min] + cvse[min]
  within1se <- which(cve >= l.se & cve <= u.se)
  min1se <- which.max(lambda %in% lambda[within1se])

  # Bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit,
              min=min, lambda.min=lambda[min],
              min1se = min1se, lambda.1se = lambda[min1se],
              null.dev=mean(loss.plmm(y, rep(mean(y), n))), Bias=Bias, Loss = E)
  if (returnY) val$Y <- Y
  structure(val, class="cv.plmm")
}
