#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param estimated_V Estimated variance-covariance matrix using all observations when computing BLUP; NULL if type = "lp" in cv.plmm. 
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
  cv.args$prep$U <- full_cv_prep$U[fold!=i, , drop=FALSE]
  cv.args$prep$y <- full_cv_prep$y[fold!=i] 
  
  # extract test set (comes from cv prep on full data)
  test_X <- full_cv_prep$std_X[fold==i, , drop=FALSE] 
  test_y <- full_cv_prep$y[fold==i]
  
  # OLD WAY 
  # NB: eta used in each fold comes from the overall fit.args. If user-supplied, then 
  # use that in all fold; if not, estimate eta in each fold 
  
  # NEW WAY 
  # I moved estimate_eta() into prep, so that this is only done once. In doing this,
  # I am assuming that the eta is the same across the training and testing data.

  # fit a plmm within each fold 
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  fit.i <- do.call("plmm_fit", cv.args)

  if(type == "lp"){
    yhat <- predict.list(fit = fit.i,
                          oldX = cv.args$prep$std_X,
                          newX = test_X,
                          type = 'lp')
    # yhat <- matrix(data = drop(Xbeta), nrow = length(y_test)) 
  }
  
  if (type == 'blup'){
    # estimated_V here comes from the overall fit in cv_plmm.R, an n*n matrix 
    V21 <- estimated_V[fold==i, fold!=i, drop = FALSE] 
    V11 <- estimated_V[fold!=i, fold!=i, drop = FALSE] 
    
    yhat <- predict.list(fit = fit.i,
                         oldX = cv.args$prep$std_X,
                         newX = test_X,
                         type = 'blup',
                         V11 = V11,
                         V21 = V21, ...)
    
  }

  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
