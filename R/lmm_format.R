#' LMM format: a function to format the output of a model constructed with \code{lmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{lmm_fit}
#' @param dfmax dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. 

#' 
#' @return A list with the components: 
#' * beta_vals: The estimated beta values at each value of lambda
#' * eta: The estimated eta value 
#' * ns_idx: COME BACK HERE
#' * iter: The number of iterations at each value of lambda (MAYBE take this out)
#' * converged: The convergence status at each value of lambda
#' * X: if returnX = TRUE and size SUX < 100 Mb, the original X will be returned
#' 
#' 
#' @keywords internal
#'
lmm_format <- function(fit,
                        dfmax = fit$ncol_X + 1, 
                        X){
 
  # reverse the transformations of the beta values 
  beta_vals <- lmm_untransform(res_b = fit$res$coefficients,
                           ns = fit$ns,
                           ncol_X = fit$p,
                           std_X = fit$std_X,
                           SUX = fit$rot_X,
                           std_SUX = fit$stdrot_X)
  
  if(fit$trace){cat("\nBeta values are estimated -- almost done!\n")}
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows
  names(beta_vals) <- c("(Intercept)", fit$snp_names)
  
  ## output
  val <- structure(list(beta_vals = beta_vals,
                        lm.fit = fit$res,
                        eta = fit$eta,
                        s = fit$s,
                        U = fit$U,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        p = fit$p,
                        n = fit$n,
                   class = "lmm"))
  
  return(val)
  
  
}

