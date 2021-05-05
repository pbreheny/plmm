#' Fit a linear model with lasso regularization
#'
#' This function allows you to fit a linear model via lasso-penalized maximum likelihood with the same output values as those of lmm_lasso_eta. Primarily used for simulation and comparison with lmm_lasso_eta.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param ... Additional optional arguments
#' @importFrom glmnet glmnet
#' @export

lasso <- function(X, y, p1, ...) {
  args.list <- list(...)
  args.list$x <- X
  args.list$y <- y
  fit <- do.call('glmnet', args.list)
  sel <- sapply(stats::predict(fit, type='nonzero'), length)
  if (!(p1 %in% sel)){
    args.list$nlambda <- length(fit$lambda)
    iter <- 1
    while (!(p1 %in% sel) & iter < 10){
      nlambda <- args.list$nlambda
      args.list$nlambda <- (nlambda + 10)
      fit <- do.call('glmnet', args.list)
      sel <- sapply(stats::predict(fit, type='nonzero'), length)
      iter <- iter + 1
    }
  }
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1, 1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = NA,
              eta = NA))
}
