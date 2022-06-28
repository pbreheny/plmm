#' Fit a penalizedLMM with lasso regularization
#'
#' This function allows you to fit a linear mixed model via lasso-penalized maximum likelihood.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param K Matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects. Default is \eqn{\frac{1}{p} XX^T}
#' @importFrom zeallot %<-%
#' @export



lmm_lasso <- function(X, y, p1, standardize = FALSE, K = NULL) {
  S <- U <- eta <- NULL
  if (is.null(K)){
    c(S, U, eta) %<-% plmm_null(X, y)
  } else {
    c(S, U, eta) %<-% plmm_null(K, y)
  }
  W <- diag((eta * S + (1 - eta))^(-1/2))
  SUX <- W %*% crossprod(U, cbind(1, X))
  SUy <- drop(W %*% crossprod(U, y))
  fit <- glmnet::glmnet(SUX, SUy, standardize = standardize, intercept = FALSE, penalty.factor = c(0, rep(1, ncol(X))))
  sel <- sapply(stats::predict(fit, type='nonzero'), length)
  coef <- coef(fit, min(fit$lambda[sel <= (p1 + 1)]))[-c(1:2)] # this is to remove the empty intercept, and the first coefficient, which corresponds to the rotated manual intercept
  names(coef) <- colnames(X)
  coef_pred <- coef(fit, min(fit$lambda[sel <= (p1 + 1)]))[-1] # this is to remove the empty intercept - return rotated int
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              coef_pred = coef_pred,
              delta = (1/eta) - 1,
              eta = eta))
}
