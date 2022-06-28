#' Fit a penalizedLMM with lasso regularization using penalizedLMM
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param K Similarity matrix. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param ... Additional optional arguments
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples
#' RRM <- genRelatednessMat(X = scale(admix$X))
#' fit <- plmm_lasso(X = admix$X, y = admix$y, K = RRM, p1 = 10)

plmm_lasso <- function(X, y, K, p1, ...) {
  args.list <- list(...)
  args.list$X <- X
  args.list$y <- y
  args.list$K <- K
  args.list$penalty <- 'lasso'
  fit <- do.call('plmm', args.list)
  sel <- predict.plmm(fit, type = "nvar")
  if (!(p1 %in% sel)){
    args.list$lambda.min <- ifelse(nrow(X) > ncol(X), 0.001, 0.05)
    lams <- seq(from = (args.list$lambda.min + 1e-4), to = .999, length.out = 100)
    iter <- 1
    while (!(p1 %in% sel) & iter <= 100){
      args.list$lambda.min <- lams[iter]
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




