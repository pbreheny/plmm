#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#' 
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param X Design matrix. May include clinical covariates and other non-SNP data. 
#' 
#' @returns A list with the components: 
#'  * `beta_vals`: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * `rotated_scale_beta_vals`: the matrix of estimated coefficients on the ~rotated~ scale. This is the scale on which the model was fit. 
#'  * `lambda`: a numeric vector of the lasso tuning parameter values used in model fitting. 
#'  * `eta`: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
#'  * `s`: a vectof of the eigenvalues of relatedness matrix `K`; see `relatedness_mat()` for details.
#'  * `U`: a matrix of the eigenvalues of relatedness matrix `K`
#'  * `rot_y`: the vector of outcome values on the rotated scale. This is the scale on which the model was fit. 
#'  * `linear.predictors`: the matrix resulting from the product of `stdrot_X` and the estimated coefficients on the ~rotated~ scale.
#'  * `penalty`: character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * `gamma`: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
#'  * `alpha`: numeric value indicating the elastic net tuning parameter. 
#'  * `convex.min`: NULL (This is an option we will add in the future!)
#'  * `loss`: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `penalty.factor`: vector of indicators corresponding to each predictor, where 1 = predictor was penalized. 
#'  * `ns_idx`: vector with the indicies of predictors which were constant features (i.e., had no variation).
#'  * `p`: the number of features 
#'  * `n`: the number of observations (instances)
#'  * `iter`: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: vector of logical values indicating whether the model fitting converged at each value of `lambda`
#' 
#' @keywords internal
#'

plmm_format <- function(fit, X){
  
  # reverse the transformations of the beta values 
  if(fit$trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  
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
    
  
