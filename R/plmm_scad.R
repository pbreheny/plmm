#' Fit a penalizedLMM with SCAD regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param K Similarity matrix. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' RRM <- relatedness_mat(X = scale(admix$X))
#' fit <- plmm_scad(X = admix$X, y = admix$y, p1 = 5, K = RRM) #FIXME: throws Error in seq.default(log(lambda.max), log(lambda.min * lambda.max), len = nlambda) : 'from' must be a finite number



plmm_scad <- function(X, y, p1, standardize = FALSE, K) {
  if (missing(K)){ # case 1: use default K 
    fit <- plmm(X, y, penalty = 'SCAD', standardizeX = standardize)
  } else { # case 2: use user-supplied K 
    fit <- plmm(X, y, K, penalty = 'SCAD', standardizeX = standardize)
  }
  sel <- predict.plmm(fit, type = "nvar")
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = (1/fit$eta) - 1,
              eta = fit$eta))
}


