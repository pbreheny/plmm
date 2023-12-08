#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
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

plmm_format <- function(prep, fit, X){
  
  # reverse the transformations of the beta values 
  if(prep$trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  
  # get beta values back in original scale; reverse the PRE-ROTATION standardization 
  untransformed_b2 <- untransform(untransformed_b1 = fit$untransformed_b1,
                                  ns = fit$ns,
                                  std_X_details = fit$std_X_details,
                                  p = fit$p)
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows, lambda values on the columns
  varnames <- c("(Intercept)", fit$snp_names)
  dimnames(untransformed_b2) <- list(varnames, lamNames(fit$lambda))
  
  colnames(fit$linear.predictors) <- lamNames(fit$lambda)

  ## output
  val <- structure(list(beta_vals = untransformed_b2,
                        rotated_scale_beta_vals = fit$untransformed_b1,
                        lambda = fit$lambda,
                        eta = fit$eta,
                        s = fit$s,
                        U = fit$U,
                        rot_y = fit$rot_y,
                        linear.predictors = fit$linear.predictors,
                        penalty = fit$penalty,
                        gamma = fit$gamma,
                        alpha = fit$alpha,
                        convex.min = fit$convex.min,
                        loss = fit$loss,
                        penalty.factor = fit$penalty.factor,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        p = fit$p,
                        n = fit$n, 
                        iter = fit$iter,
                        converged = fit$converged),
                   class = "plmm")
  
  
  return(val)
  
  
}
    
  
