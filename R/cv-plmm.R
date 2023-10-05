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
                      trace = trace,
                      ...)) # ... additional arguments to plmm_prep()

  prep <- do.call('plmm_prep', prep.args)
  
  
  # implement full model fit 
  fit.args <- c(list(prep = prep, penalty = penalty), list(...))
  fit <- do.call('plmm_fit', fit.args)
  fit_to_return <- plmm_format(fit, X = X)
  
  # set up arguments for cv 
  cv.args <- fit.args
  cv.args$warn <- FALSE
  cv.args$lambda <- fit$lambda 

  # initialize objects to hold CV results 
  n <- length(fit$y)
  cve_mat <- cvse_mat <- matrix(NA, nrow=nfolds, ncol=length(fit$lambda))
  E <- Y <- array(dim = c(length(fit$y), length(fit$lambda), nfolds))
  
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
          ")
    }
    fold <- sample(1:nrow(prep$stdrot_X) %% nfolds)
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
  # NB: I am splitting the rows of stdrot_X into folds, so I am indexing 'k', not 'n'
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
                 K = K)
      if (trace) {setTxtProgressBar(pb, i)}

    }
    
    # fill in array with loss and yhat values 
    E[,,i] <- res$loss
    Y[,,i] <- res$yhat
    
    
    # eliminate saturated lambda values, if any
    ind <- which(apply(is.finite(res$loss), 2, all)) # index for lambda values to keep
    res$loss <- res$loss[, ind, drop=FALSE]
    res$yhat <- res$yhat[, ind]
    lambda <- fit$lambda[ind]
    
    # assess loss 
    cve <- apply(res$loss, 2, mean)
    cvse <- apply(res$loss, 2, stats::sd) / sqrt(length(fit$y))
    
    cve_mat[i, 1:res$nl] <- cve
    cvse_mat[i, 1:res$nl] <- cvse
  }

  # return min lambda idx
  avg_cve_over_folds <- colSums(cve_mat)/nrow(cve_mat)
  min <- which.min(avg_cve_over_folds)

  # return lambda 1se idx
  avg_cvse_over_folds <- colSums(cvse_mat)/nrow(cvse_mat)
  l.se <- avg_cve_over_folds[min] - avg_cvse_over_folds[min]
  u.se <- avg_cve_over_folds[min] + avg_cvse_over_folds[min]
  within1se <- which(avg_cve_over_folds >= l.se & avg_cve_over_folds <= u.se)
  min1se <- which.max(lambda %in% lambda[within1se])

  # bias correction
  if(returnBiasDetails){
    Bias <- mean(cve_mat[,min] - apply(cve_mat, 2, min))
  }
  

  val <- list(type=type, 
              avg_cve_over_folds=avg_cve_over_folds,
              avg_cvse_over_folds=avg_cvse_over_folds,
              cve_mat = cve_mat,
              cvse_mat = cvse_mat,
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
