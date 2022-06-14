#' Fit a penalizedLMM with MCP regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param X_for_K X matrix used to compute the similarity matrix, K. Defaults to XX^T. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @importFrom zeallot %<-%
#' @export



plmm_mcp <- function(X, y, p1, standardize = FALSE, X_for_K = NULL) {
  if (is.null(X_for_K)){
    fit <- plmm(X = X,
                y = y,
                V = X%*%t(X), # TODO: verify this is appropriate (changed on June 14, 2022)
                penalty = 'MCP',
                standardizeX = standardize)
  } else {
    fit <- plmm(X = X,
                y = y,
                V = X_for_K,
                penalty = 'MCP',
                standardizeX = standardize)
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


