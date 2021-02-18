#' Cross-validation internal function for cv.plmm that computes prediction based on the unrotated data.
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param XX Transformed design matrix including intercept. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied with a  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param X_unrotated Non-transformed design matrix, which should not include an intercept. May include clinical covariates and other non-SNP data.
#' @param y_unrotated Non-transformed outcome vector.
#' @param fold n-length vector of fold-assignments.
#' @param cv.args List of additional arguments to be passed to plmm.
#' @export

cvf_predict_unrotated <- function(i, XX, y, X_unrotated, y_unrotated, fold, cv.args) {
  # fit on the rotated data
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("plmm", cv.args)

  # predict using the non-rotated data
  X2 <- cbind(1, X_unrotated[fold==i, , drop=FALSE])
  y2 <- y_unrotated[fold==i]
  yhat <- matrix(predict.plmm(fit.i, X=X2, type="response", lambda=fit.i$lambda), length(y2))
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
