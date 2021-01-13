#' Fit a linear model with lasso regularization
#'
#' This function allows you to fit a linear model via lasso-penalized maximum likelihood with the same output values as those of lmm_lasso_eta. Primarily used for simulation and comparison with lmm_lasso_eta.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @export

lasso <- function(X, y, p1, standardize = FALSE) {
  fit <- glmnet::glmnet(X, y, standardize = standardize)
  sel <- sapply(stats::predict(fit, type='nonzero'), length)
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[,1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)) - 1,
              coef = coef,
              delta = NA,
              eta = NA))
}
