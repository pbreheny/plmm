#' Fit a penalizedLMM with MCP regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param K Similarity matrix. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param ... Other arguments to plmm()
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' plmm_mcp_fit <- plmm_mcp(admix$X, admix$y, p1 = 5)
#'

plmm_mcp <- function(X, y, p1, standardize = FALSE, K, ...) {
  if (missing(K)){ # case 1: no K supplied, so use the plmm() default: 
    fit <- plmm(X = X,
                y = y, 
                penalty = 'MCP',
                ...)
  } else { # case 2: user-specified K
    fit <- plmm(X = X,
                y = y,
                K = K,
                penalty = 'MCP',
                ...)
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


