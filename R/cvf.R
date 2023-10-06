#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param ... Optional arguments to `predict.list`
#' @importFrom zeallot %<-%
#' @keywords internal
cvf <- function(i, fold, type, cv.args, K, ...) {

  # save the 'prep' object from the plmm_prep() in cv.plmm
  full_cv_prep <- cv.args$prep
  
  # subset std_X, and y to match fold indices 
  #   (and in so doing, leave out the ith fold)
  cv.args$prep$rot_y <- full_cv_prep$rot_y[fold!=i] 
  cv.args$prep$stdrot_X <- full_cv_prep$stdrot_X[fold!=i, ,drop=FALSE]
  # attr(cv.args$prep$stdrot_X, 'scale') <- attr(full_cv_prep$stdrot_X, 'scale')[fold!=i]
  # TODO: do I need center & scale values above? What to pass to untransform() via predict.list()?
  
  # for estimating V (v_hat()) in BLUP case, I will need to get U and s for the ith fold:
  cv.args$prep$U <- full_cv_prep$U[,fold!=i,drop=FALSE] # TODO: generalize this to accommodate p > n case 
  cv.args$prep$s <- full_cv_prep$s[fold!=i]
  
  # extract test set (comes from cv prep on full data)
  # NB: test data is on *standardized* scale 
  # test_U <- full_cv_prep$U[,fold==i,drop=FALSE]
  # test_Vt <- full_cv_prep$Vt[,fold==i,drop=FALSE]
  # test_d <- sqrt(full_cv_prep$s[fold==i]*full_cv_prep$p)
  # test_X <- test_U%*%diag(test_d)%*%test_Vt
  test_X <- full_cv_prep$std_X
  test_y <- full_cv_prep$y # TODO: think more about this choice... 
  
  # NB: eta used in each fold comes from the overall fit.args. If user-supplied, then 
  # use that in all fold; if not, estimate eta in each fold 

  # fit a plmm within each fold 
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  fit.i <- do.call("plmm_fit", cv.args)
  
  if(type == "lp"){
    yhat <- predict.list(fit = fit.i,
                          std_X = cv.args$prep$std_X,
                          newX = test_X,
                          type = 'lp',
                         stdrot_X_scale = attr(full_cv_prep$stdrot_X, 'scale'))
    # yhat <- matrix(data = drop(Xbeta), nrow = length(y_test)) 
  }
  
  if (type == 'blup'){

    V11 <- v_hat(fit.i, K = K)
    V21 <- v_hat(fit.i, K = K)
    
    yhat <- predict.list(fit = fit.i,
                         std_X = cv.args$prep$std_X,
                         newX = test_X,
                         type = 'blup',
                         stdrot_X_scale = attr(full_cv_prep$stdrot_X, 'scale'),
                         prep = cv.args$prep, 
                         V11 = V11,
                         V21 = V21, ...)
    
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
