#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'response'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv.args List of additional arguments to be passed to plmm.
#' @param ... Optional arguments to `predict.list`
#' @importFrom zeallot %<-%
#' @keywords internal
cvf <- function(i, fold, type, cv.args, ...) {
 
  # subset std_X, and y to match fold indices 
  #   (and in so doing, leave out the ith fold)
  X_train <- cv.args$prep$std_X[fold!=i, ,drop=FALSE] 
  y_train <- cv.args$prep$y[fold!=i] 
  K_train <- cv.args$K[fold!=i, fold!=i, drop=FALSE]
  
  # V <- cv.args$estimated_V
  # V_train <- V[fold!=i, fold!=i]
  # V_train_test <- V[fold!=i, fold=i, drop = FALSE]   
  
  # do SVD inside each fold using training data 
  prep.args.i <- c(list(X = X_train,
                        y = y_train,
                        K = K_train,
                        returnX = cv.args$prep$returnX,
                        trace = cv.args$prep$trace))

  
  prep.i <- do.call('plmm_prep', prep.args.i)
  
  # fit a plmm within each fold 
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm.R line 63 
  fit.args.i <- c(list(prep.i, penalty = cv.args$penalty, lambda = cv.args$lambda), list(...))
  fit.i <- do.call("plmm_fit", fit.args.i)
  ns.i <- attr(fit.i$std_X, "nonsingular")
  
  # extract test set
  X_test <- cv.args$prep$std_X[fold==i, ns.i, drop=FALSE] 
  # NOTE: make sure to align columns of X_test with the NONSINGULAR columns of std(X_train)
  y_test <- cv.args$prep$y[fold==i]
  K_test <- cv.args$K[fold==i, fold==i, drop=FALSE]
  
  
  beta <- coef.list(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- predict.list(fit = fit.i, newX = X_test, type = 'response', lambda = fit.i$lambda)
  yhat <- matrix(data = drop(Xbeta), nrow = length(y_test))
  
  
  if (type == 'blup'){
    yhat <- predict.list(fit = fit.i, newX = X_test, type = 'blup',
                         lambda = fit.i$lambda, prep = cv.args$prep, ...)
    
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y_test, yhat[,ll]))
  # bias <- cv_bias(X_train, X_test, y_train, y_test, V_train, V_train_test)
  list(loss=loss,
       # bias=bias,
       nl=length(fit.i$lambda),
       yhat=yhat)
}