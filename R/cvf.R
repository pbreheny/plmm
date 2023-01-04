#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param XX Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied with a matrix containing only SNP data for computation of GRM.
#' @param y Original continuous outcome vector.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'response'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv.args List of additional arguments to be passed to plmm.
#' @param ... Optional arguments
#' @importFrom zeallot %<-%
#' @export
#' 



cvf <- function(i, XX, y, fold, type, cv.args, ...) {
  # overwrite X and y arguments to leave out the ith fold
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  # NB: inside each fold, we are not re-doing the prep steps like SVD, rotation, etc.
  fit.i <- do.call("plmm_fit", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i]

  beta <- coef.plmm(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- predict.plmm(fit.i, newX = X2, type = 'response', lambda = fit.i$lambda)
  yhat <- matrix(drop(Xbeta), length(y2))
  
  if (type == 'blup'){
    yhat <- predict.plmm(fit.i, newX = X2, type = 'blup', lambda = fit.i$lambda,
                         X = cv.args$X, y = cv.args$y, U = fit.i$U, S = fit.i$S,
                         eta = fit.i$eta, 
                         covariance = K[fold == i, fold != i, drop = FALSE])
  }
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
