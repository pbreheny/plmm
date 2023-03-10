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
cvf <- function(i, fold, type, cv.args, estimated_V, ...) {
  
  # save the 'prep' object from the plmm_prep() in cv.plmm
  full_cv_prep <- cv.args$prep
  
  # subset std_X, and y to match fold indices 
  #   (and in so doing, leave out the ith fold)
  cv.args$prep$std_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]
  cv.args$prep$U <- full_cv_prep$U[fold!=i, ,drop=FALSE]
  cv.args$prep$y <- full_cv_prep$y[fold!=i] 
  
  # eta used in each fold comes from the overall fit.args. If user-supplied, then use that in all fold; if not, estimate eta in each fold 
  
  # create copies of the training sets (to be used in cv_bias)
  # TODO: could this be more efficient/readable? 
  X_train <- full_cv_prep$std_X[fold!=i, ,drop=FALSE] 
  y_train <- full_cv_prep$y[fold!=i] 
  
  # extract test set (comes from cv prep on full data)
  X_test <- full_cv_prep$std_X[fold==i, , drop=FALSE] 
  y_test <- full_cv_prep$y[fold==i]
  
  # fit a plmm within each fold 
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm.R line 63 
  fit.i <- do.call("plmm_fit", cv.args)
  
  beta <- coef.list(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- predict.list(fit = fit.i, newX = X_test, type = 'response', lambda = fit.i$lambda)
  yhat <- matrix(data = drop(Xbeta), nrow = length(y_test))
  
  if (type == 'blup'){
    # used as V21 when computing blup (V_test_train)
    V21 <- estimated_V[fold==i, fold!=i, drop = FALSE] 
    V11 <- estimated_V[fold!=i, fold!=i, drop = FALSE] 
    
    yhat <- predict.list(fit = fit.i, newX = X_test, type = 'blup',
                         lambda = fit.i$lambda, prep = cv.args$prep, 
                         V11 = V11, V21 = V21, ...)
    
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y_test, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
