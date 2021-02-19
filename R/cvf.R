#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param XX Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param fold n-length vector of fold-assignments.
#' @param cv.args List of additional arguments to be passed to plmm.
#' @param predict.args List of arguments required for prediction. Varies according to prediction type (response vs. individual).
#' @export

cvf <- function(i, XX, y, fold, cv.args, predict.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("plmm", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  predict.args$object <- fit.i
  predict.args$newX <- X2
  if (predict.args$type == 'individual'){
    predict.args$X <- cv.args$X
    predict.args$y <- cv.args$y
    predict.args$eigenvectors <- fit.i$eigenvectors
    predict.args$eigenvalues <- fit.i$eigenvalues
    predict.args$eta <- fit.i$eta
    predict.args$covariance <- tcrossprod(ncvreg::std(X2), ncvreg::std(cv.args$X))/ncol(X2)
  }
  # yhat <- matrix(predict.plmm(fit.i, newX=X2, type="response", lambda=fit.i$lambda), length(y2))
  yhat <- matrix(do.call("predict.plmm", predict.args), length(y2))
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  # if (max(loss) > 1e10) browser()
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
