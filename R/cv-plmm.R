#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized penalized linear mixed models over a grid of values for the regularization parameter lambda.
#' @param X Design matrix for model fitting. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector for model fitting.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix, in the form of (1) the relatedness matrix estimated from the data (default), (2) a user-supplied matrix, or (3) a user-supplied list with components 'd' and 'u' as created by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#' @param eta_star Optional arg. to \code{plmm_prep}. Defaults to NULL.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param penalty.factor Optional arg. to \code{plmm_prep}. Defaults to 1 for all predictors (except the intercept). 
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'blup'} predictions are based on the *sum* of the linear predictor and the estimated random effect (BLUP). Defaults to 'lp'.
#' @param ... Additional arguments to \code{plmm_fit}
#' @param cluster cv.plmm can be run in parallel across a cluster using the parallel package. The cluster must be set up in advance using the makeCluster function from that package. The cluster must then be passed to cv.plmm.
#' @param nfolds The number of cross-validation folds. Default is 10.
#' @param fold Which fold each observation belongs to. By default the observations are randomly assigned.
#' @param seed You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY Should cv.plmm return the linear predictors from the cross-validation folds? Default is FALSE; if TRUE, this will return a matrix in which the element for row i, column j is the fitted value for observation i from the fold in which observation i was excluded from the fit, at the jth value of lambda.
#' @param returnBiasDetails Logical: should the cross-validation bias (numeric value) and loss (n x p matrix) be returned? Defaults to FALSE. 
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @export
#' 
#' @examples 
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, seed = 321)
#' \dontrun{
#' cv_s <- summary.cv.plmm(cv_fit, lambda = "1se")
#' print(cv_s)
#' plot(cv_fit)
#' }
#' 
#' 

cv.plmm <- function(X,
                    y,
                    k = NULL,
                    K = NULL,
                    diag_K = NULL,
                    eta_star = NULL,
                    penalty = c("MCP", "SCAD", "lasso"),
                    penalty.factor = rep(1, ncol(X)),
                    type = 'lp',
                    ...,
                    cluster,
                    nfolds=10,
                    seed,
                    fold = NULL,
                    returnY=FALSE,
                    returnBiasDetails = FALSE,
                    trace=FALSE) {


  # default type is 'lp'
  if(missing(type)) {type == 'lp'} 

  # determine penalty 
  penalty <- match.arg(penalty)
  
  # implement preparation steps for model fitting 
  prep.args <- c(list(X = X,
                      y = y,
                      k = k,
                      K = K,
                      diag_K = diag_K,
                      eta_star = eta_star,
                      penalty.factor = penalty.factor,
                      trace,
                      ...)) # ... additional arguments to plmm_prep()
  
  prep <- do.call('plmm_prep', prep.args)
  
  
  # implement full model fit 
  fit.args <- c(list(prep = prep, penalty = penalty), list(...))
  fit <- do.call('plmm_fit', fit.args)
  fit_to_return <- plmm_format(fit, X = X, K = K)
  
  # set up arguments for cv 
  cv.args <- fit.args
  cv.args$warn <- FALSE
  cv.args$lambda <- fit$lambda 
  
  estimated_V <- NULL 
  if (type == 'blup') {
    estimated_V <- fit$eta * tcrossprod(fit$U %*% diag(fit$S), fit$U) + (1-fit$eta)*diag(nrow = nrow(fit$U)) 
  }

  # initialize objects to hold CV results 
  n <- length(fit$y)
  E <- Y <- matrix(NA, nrow=nrow(X), ncol=length(fit$lambda))

  
  # set up folds for cross validation 
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  
  sde <- sqrt(.Machine$double.eps)
  
  if(is.null(fold)) {
    if(trace){
      cat("'Fold' argument is either NULL or missing; assigning folds randomly (by default). 
          \nTo specify folds for each observation, supply a vector with fold assignments.")
    }
    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  
  # set up cluster if user-specified
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "K", "fold", "type", "cv.args", "estimated_V"), envir=environment())
    parallel::clusterCall(cluster, function() library(penalizedLMM))
    fold.results <- parallel::parLapply(cl=cluster, X=1:max(fold), fun=cvf, X=X, y=y, fold=fold, type=type, cv.args=cv.args, 
                                        estimated_V = estimated_V)
  }
  
  if (trace) cat("\nStarting cross validation\n")  
  # set up progress bar -- this can take a while
  if(trace){pb <- txtProgressBar(min = 0, max = nfolds, style = 3)}
  for (i in 1:nfolds) {
    # case 1: user-specified cluster
    if (!missing(cluster)) {
      res <- fold.results[[i]] # refers to lines from above
      if (trace) {setTxtProgressBar(pb, i)}
    } else {
      # case 2: cluster NOT user specified 
      res <- cvf(i = i,
                 fold = fold,
                 type = type,
                 cv.args = cv.args,
                 estimated_V = estimated_V)
      if (trace) {setTxtProgressBar(pb, i)}
    }
    # update E and Y
    E[fold==i, 1:res$nl] <- res$loss
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  # eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all)) # index for lambda values to keep
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  # return min lambda idx
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(n)
  min <- which.min(cve)

  # return lambda 1se idx
  l.se <- cve[min] - cvse[min]
  u.se <- cve[min] + cvse[min]
  within1se <- which(cve >= l.se & cve <= u.se)
  min1se <- which.max(lambda %in% lambda[within1se])

  # bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(type=type, 
              cve=cve,
              cvse=cvse,
              fold=fold,
              lambda=lambda,
              fit=fit_to_return,
              min=min,
              lambda.min=lambda[min],
              min1se = min1se,
              lambda.1se = lambda[min1se],
              null.dev=mean(loss.plmm(y, rep(mean(y), n))))
  if (returnY) val$Y <- Y
  if (returnBiasDetails){
    val$Bias <- Bias
    val$Loss <- E
  }
  structure(val, class="cv.plmm")
}
