#' Fit a penalizedLMM with lasso regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param ... Additional optional arguments
#' @importFrom zeallot %<-%
#' @export



plmm_lasso <- function(X, y, p1, ...) {
  args.list <- list(...)
  args.list$X <- X
  args.list$y <- y
  args.list$penalty <- 'lasso'
  fit <- do.call('plmm', args.list)
  sel <- predict.plmm(fit, type = "nvar")
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  coef_pred <- coef(fit, min(fit$lambda[sel <= p1]))
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              coef_pred = coef_pred,
              delta = (1/fit$eta) - 1,
              eta = fit$eta))
}


