#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param convex convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. 
#' 
#' @return A list with the components: 
#' * beta_vals: The estimated beta values at each value of lambda
#' * eta: The estimated eta value 
#' * lambda: The sequence of lambda values used in model fitting
#' * penalty: A string indicating the type of penalty used to fit the model
#' * ns_idx: COME BACK HERE
#' * iter: The number of iterations at each value of lambda (MAYBE take this out)
#' * converged: The convergence status at each value of lambda
#' 
#' 
#' @keywords internal
#'
#' @examples
#' 
#' \dontrun{
#' # NB: this is an internal function (i.e. it's not exported)
#' # these examples won't work unless you have plmm::: 
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
#' fit1 <- plmm_fit(prep = prep1)
#' res1 <- plmm_format(fit1)
#' 
#' }

plmm_format <- function(fit,
                        convex = TRUE,
                        dfmax = fit$p + 1, 
                        X){
  # eliminate saturated lambda values, if any
  ind <- !is.na(fit$iter)
  iter <- fit$iter[ind]
  converged <- fit$converged[ind]
  lambda <- fit$lambda[ind]
  loss <- fit$loss[ind]
  if (fit$warn & sum(iter) == fit$max.iter) warning("\nMaximum number of iterations reached")
  convex.min <- if (convex) convexMin(b = fit$b,
                                      X = fit$stdrot_X,
                                      penalty = fit$penalty,
                                      gamma = fit$gamma, 
                                      l2 = fit$lambda*(1-fit$alpha),
                                      family = 'gaussian',
                                      penalty.factor = fit$penalty.factor) else NULL
  
  # reverse the transformations of the beta values 
  beta_vals <- untransform(res_b = fit$b,
                           ns = fit$ns,
                           p = fit$p,
                           std_X = ncvreg::std(X),
                           rot_X = fit$rot_X,
                           stdrot_X = fit$stdrot_X)
  
  if(fit$trace){cat("\nBeta values are estimated -- almost done!")}
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows, lambda values on the columns
  varnames <- c("(Intercept)", fit$snp_names)
  dimnames(beta_vals) <- list(varnames, lamNames(fit$lambda))
  
  colnames(fit$linear.predictors) <- lamNames(fit$lambda)

  ## output
  val <- structure(list(beta_vals = beta_vals,
                        lambda = lambda,
                        eta = fit$eta,
                        s = fit$s,
                        U = fit$U,
                        rot_y = fit$rot_y,
                        linear.predictors = fit$linear.predictors,
                        penalty = fit$penalty,
                        gamma = fit$gamma,
                        alpha = fit$alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = fit$penalty.factor,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        p = fit$p,
                        n = fit$n, 
                        iter = iter,
                        converged = converged),
                   class = "plmm")
  
  
  return(val)
  
  
}
    
  
