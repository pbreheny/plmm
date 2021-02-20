#' Cross-validation internal function for cv.plmm
#'
#' Internal function for cv.plmm which calls plmm on a fold subset of the original data.
#' @param i Fold number to be excluded from fit.
#' @param XX Design matrix. May include clinical covariates and other non-SNP data. If this is the case, XX_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param XX_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param y Original continuous outcome vector.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'response'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv.args List of additional arguments to be passed to plmm.
#' @export

cvf <- function(i, XX, XX_for_K=XX, y, fold, type, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$X_for_K <- XX_for_K[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("plmm", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  X_for_K2 <- XX_for_K[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  if (type == 'response'){
    yhat <- matrix(predict.plmm(object=fit.i, newX=X2, type=type, lambda=fit.i$lambda), length(y2))
  } else if (type == 'individual'){
    covariance <- tcrossprod(ncvreg::std(X_for_K2), ncvreg::std(cv.args$X_for_K))/ncol(X_for_K2)
    yhat <- matrix(predict.plmm(object=fit.i, newX=X2, type=type, lambda=fit.i$lambda,
                                XX=cv.args$X, y=cv.args$y, U=fit.i$U, S=fit.i$S,
                                eta=fit.i$eta, covariance=covariance),
                   length(y2))
  }

  # predict.args <- list()
  # predict.args$type <- type
  # predict.args$lambda <- fit.i$lambda
  # predict.args$fit <- fit.i
  # predict.args$newX <- X2
  # if (predict.args$type == 'individual'){
  #   predict.args$X <- cv.args$X
  #   predict.args$y <- cv.args$y
  #   predict.args$U <- fit.i$U
  #   predict.args$S <- fit.i$S
  #   predict.args$eta <- fit.i$eta
  #   predict.args$covariance <- tcrossprod(ncvreg::std(X2), ncvreg::std(cv.args$X))/ncol(X2)
  # }
  # yhat <- matrix(do.call("predict.plmm", predict.args), length(y2))
  loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  # if (max(loss) > 1e10) browser()
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
