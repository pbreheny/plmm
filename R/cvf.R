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
  
  # save the 'prep' object from the plmm_prep() in cv.plmm
  full_cv_prep <- cv.args$prep
  
  # subset std_X, and y to match fold indices 
  #   (and in so doing, leave out the ith fold)
  cv.args$prep$std_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]
  cv.args$prep$U <- full_cv_prep$U[fold!=i, ,drop=FALSE]
  cv.args$prep$y <- full_cv_prep$y[fold!=i]
  cv.args$eta <- NULL # we will re-estimate eta in each fold
  cv.args$estimated_V <- NULL # TODO: revisit this later 
  
  # create copies of the training sets (to be used in cv_bias)
  # TODO: could this be more efficient/readable? 
  X_train <- full_cv_prep$std_X[fold!=i, ,drop=FALSE] 
  y_train <- full_cv_prep$y[fold!=i] 
  
  # extract test set (comes from cv prep on full data)
  X_test <- full_cv_prep$std_X[fold==i, , drop=FALSE] 
  y_test <- full_cv_prep$y[fold==i]
  
  
  # V <- cv.args$estimated_V
  # V_train <- V[fold!=i, fold!=i]
  # V_train_test <- V[fold!=i, fold=i, drop = FALSE]   
  
  # fit a plmm within each fold 
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm.R line 63 
  fit.i <- do.call("plmm_fit", cv.args)
  
  beta <- coef.list(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- predict.list(fit = fit.i, newX = X_test, type = 'response', lambda = fit.i$lambda)
  yhat <- matrix(data = drop(Xbeta), nrow = length(y_test))
  
  if (type == 'blup'){
    yhat <- predict.list(fit = fit.i, newX = X_test, type = 'blup',
                         lambda = fit.i$lambda, prep = cv.args$prep, ...)
    
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y_test, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
