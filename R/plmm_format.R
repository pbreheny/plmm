#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#' 
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param snp_names A vector of names for features, passed internally if such names are included with the data `X` passed to `plmm()`
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

plmm_format <- function(fit, snp_names = NULL){

  # reverse the transformations of the beta values ----------------------
  if (fit$trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  
  # get beta values back in original scale; reverse the PRE-ROTATION standardization 
  untransformed_b2 <- untransform(untransformed_b1 = fit$untransformed_b1,
                                  # ns = fit$ns,
                                  std_X_details = fit$std_X_details,
                                  p = fit$p)
  
  # give the matrix of beta_values readable names ----------------------------
  # SNPs (or covariates) on the rows, lambda values on the columns
  
  if (is.null(snp_names)){
    snp_names <- paste("snp", 1:fit$p, sep="") 
  } 
  
  varnames <- c("(Intercept)", snp_names) # add intercept 
  dimnames(untransformed_b2) <- list(varnames, lamNames(fit$lambda))
  colnames(fit$linear.predictors) <- lamNames(fit$lambda)
  
  # output ------------------------------------------
  val <- structure(list(beta_vals = untransformed_b2,
                        rotated_scale_beta_vals = fit$untransformed_b1,
                        lambda = fit$lambda,
                        eta = fit$eta,
                        s = fit$s,
                        U = fit$U,
                        rot_y = fit$rot_y,
                        # rot_X = fit$rot_X, # TODO: think of whether this is important to keep, e.g., for GIC
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


