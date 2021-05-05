#' Fit a linear model with lasso regularization and PC adjutment
#'
#' This function allows you to fit a linear model via lasso-penalized maximum likelihood, adjusted for the first 10 principal components (unpenalized), with the same output values as those of lmm_lasso_eta. Primarily used for simulation and comparison with lmm_lasso_eta.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param ... Additional optional arguments
#' @importFrom zeallot %<-%
#' @importFrom glmnet glmnet
#' @export


lasso_pca10 <- function(X, y, p1, ...){
  S <- U <- NULL
  K <- tcrossprod(ncvreg::std(X))/ncol(X)
  c(S, U) %<-% methods::as(eigen(K), "list")
  pc <- U[, 1:10]
  XX <- cbind(pc, X)
  args.list <- list(...)
  args.list$x <- XX
  args.list$y <- y
  args.list$penalty.factor <- rep(c(0, 1), times=c(10, ncol(X)))
  fit <- do.call('glmnet', args.list)
  sel <- sapply(stats::predict(fit, type='nonzero'), length) - 10 # remove unpenalized effects
  if (!(p1 %in% sel)){
    args.list$nlambda <- length(fit$lambda)
    iter <- 1
    while (!(p1 %in% sel)){
      nlambda <- args.list$nlambda
      args.list$nlambda <- (args.list$nlambda + 10)
      fit <- do.call('glmnet', args.list)
      sel <- sapply(stats::predict(fit, type='nonzero'), length) - 10
      iter <- iter + 1
    }
  }
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-c(1:11),1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = NA,
              eta = NA))
}
