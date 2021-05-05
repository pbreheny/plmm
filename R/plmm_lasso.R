#' Fit a penalizedLMM with lasso regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param V RRM to use
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param ... Additional optional arguments
#' @importFrom zeallot %<-%
#' @export



plmm_lasso <- function(X, y, V, p1, ...) {
  args.list <- list(...)
  args.list$X <- X
  args.list$y <- y
  args.list$V <- V
  args.list$penalty <- 'lasso'
  fit <- do.call('plmm', args.list)
  sel <- predict.plmm(fit, type = "nvar")
  if (!(p1 %in% sel)){
    args.list$nlambda <- length(fit$lambda)
    iter <- 1
    while (!(p1 %in% sel)){
      nlambda <- args.list$nlambda
      args.list$nlambda <- (nlambda + 10)
      fit <- do.call('plmm', args.list)
      sel <- predict.plmm(fit, type = "nvar")
      iter <- iter + 1
    }
  }
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = (1/fit$eta) - 1,
              eta = fit$eta))
}


