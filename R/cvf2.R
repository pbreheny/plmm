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
#' @importFrom zeallot %<-%
#' @export

cvf2 <- function(i, XX, XX_for_K=XX, y, fold, type, cv.args) {
  d <- u <- uu <- NULL
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  cv.args$X_for_K <- XX_for_K[fold!=i, , drop=FALSE]
  loss_V <- FALSE
  if(!is.null(cv.args$loss_V)){
    loss_V <- cv.args$loss_V
    cv.args$loss_V <- NULL
  }

  # if V is supplied this overrides the estimation of a similarity matrix from X_for_K
  if (!is.null(cv.args$V)){
    V <- cv.args$V
    cv.args$V <- V[fold!=i, fold!=i, drop=FALSE]
  }
  fit.i <- do.call("plmm", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i]

  beta <- coef.plmm(fit.i, fit.i$lambda, drop=FALSE) # includes intercept
  Xbeta <- cbind(1, X2) %*% beta
  yhat <- matrix(drop(Xbeta), length(y2))

  if (type == 'individual'){
    X1 <- cbind(1, cv.args$X)
    y1 <- matrix(rep(cv.args$y, ncol(beta)), nrow(X1))
    U <- fit.i$U
    S <- fit.i$S
    D_inv <- diag(nrow(U))
    if (!is.null(cv.args$V)){
      diag(D_inv) <- (S)^(-1)
      covariance <- V[fold == i, fold != i, drop = FALSE]
    } else {
      eta <- fit.i$eta
      diag(D_inv) <- (1 + eta * (S - 1))^(-1)
      X_for_K2 <- XX_for_K[fold==i, , drop=FALSE]
      covariance <- tcrossprod(ncvreg::std(X_for_K2), ncvreg::std(cv.args$X_for_K))/ncol(X_for_K2)
    }
    ranef <- covariance %*% U %*% D_inv %*% t(U) %*% (y1 - X1 %*% beta)
    yhat <- Xbeta + ranef
  }

  if (loss_V){
    c(d, u, uu) %<-% svd(V[fold == i, fold == i])
    d_inv <- diag(1/sqrt(d))
    loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(d_inv %*% t(u) %*% y2, d_inv %*% t(u) %*% yhat[,ll]))
  } else {
    loss <- sapply(1:ncol(yhat), function(ll) loss.plmm(y2, yhat[,ll]))
  }
  # if (max(loss) > 1e10) browser()
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}
