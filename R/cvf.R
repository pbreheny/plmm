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
  
  # extract test set
  X2 <- cv.args$prep$std_X[fold==i, , drop=FALSE]
  y2 <- cv.args$prep$y[fold==i]
  
  # subset std_X, U, and y to match fold indices 
  #   (and in so doing, leave out the ith fold)
  cv.args$prep$std_X <- cv.args$prep$std_X[fold!=i, ,drop=FALSE]
  cv.args$prep$U <- cv.args$prep$U[fold!=i, ,drop=FALSE]
  cv.args$prep$y <- cv.args$prep$y[fold!=i]
  cv.args$eta_star <- NULL # we will re-estimate eta in each fold
  
  # NB: inside each fold, we are not re-doing the prep steps like SVD, rotation, etc.
  fit.i <- do.call("plmm_fit", cv.args)
  
  beta <- coef.list(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- predict.list(fit = fit.i, newX = X2, type = 'response', lambda = fit.i$lambda)
  yhat <- matrix(data = drop(Xbeta), nrow = length(y2))
  
  if (type == 'blup'){
    yhat <- predict.list(fit = fit.i, newX = X2, type = 'blup',
                         lambda = fit.i$lambda, prep = cv.args$prep, ...)
                        
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
