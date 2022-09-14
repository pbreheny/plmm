#' Fit a linear mixed model with lasso regularization
#'
#' This function allows you to fit a linear mixed model via lasso-penalized maximum likelihood.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param K Matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects. Default is \eqn{\frac{1}{p} XX^T}
#' @export
#' 
#' @examples 
#' lmm_ggmix_fit <- lmm_ggmix(X = admix$X, y = admix$y, p1 = 10, K = relatedness_mat(admix$X))



lmm_ggmix <- function(X, y, p1, standardize = FALSE, K = NULL){
  if (is.null(K)){ # case 1: use default (RRM) value for K 
    fit <- ggmix::ggmix(x=X, y=y, kinship=tcrossprod(ncvreg::std(X))/ncol(X), estimation="full", standardize = standardize)
  } else { # case 2: user-specified K 
    fit <- ggmix::ggmix(x=X, y=y, kinship=K, estimation="full", standardize = standardize)
  }
  sel <- sapply(stats::predict(fit, type='nonzero'), length) - 3 # minus 3 for int, and 2 var comp.
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  coef_pred <- coef(fit, min(fit$lambda[sel <= p1]))
  eta_hat <- fit$eta[which.min(fit$lambda[sel <= p1])]
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              coef_pred = coef_pred,
              delta = (1/eta_hat) - 1,
              eta = eta_hat))
}
