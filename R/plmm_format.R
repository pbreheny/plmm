#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param convex convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' 
#' @return A list with the components: 
#' * beta_vals: The estimated beta values at each value of lambda
#' * eta: The estimated eta value 
#' * lambda: The sequence of lambda values used in model fitting
#' * penalty: A string indicating the type of penalty used to fit the model
#' * ns_idx: COME BACK HERE
#' * iter: The number of iterations at each value of lambda (MAYBE take this out)
#' * converged: The convergence status at each value of lambda
#' * X: if returnX = TRUE and size SUX < 100 Mb, the original X will be returned
#' 
#' 
#' @export
#'
#' @examples
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
#' fit1 <- plmm_fit(prep = prep1)
#' res1 <- plmm_format(fit1)

plmm_format <- function(fit,
                        convex = TRUE,
                        dfmax = p + 1){
  
  
  # eliminate saturated lambda values, if any
  ind <- !is.na(fit$iter)
  iter <- fit$iter[ind]
  converged <- fit$converged[ind]
  lambda <- fit$lambda[ind]
  loss <- fit$loss[ind]
  if (fit$warn & sum(iter) == fit$max.iter) warning("Maximum number of iterations reached")
  convex.min <- if (convex) convexMin(b = fit$b,
                                      X = fit$std_SUX,
                                      penalty = fit$penalty,
                                      gamma = fit$gamma, 
                                      l2 = fit$lambda*(1-fit$alpha),
                                      family = 'gaussian',
                                      penalty.factor = fit$penalty.factor) else NULL
  
  # reverse the transformations of the beta values 
  beta_vals <- untransform(res_b = fit$b,
                           ns = fit$ns,
                           ncol_X = fit$ncol_X,
                           std_X = fit$std_X,
                           SUX = fit$SUX,
                           std_SUX = fit$std_SUX)
  
  if(fit$trace){cat("\nBeta values are estimated -- almost done!\n")}
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows, lambda values on the columns
  varnames <- c("(Intercept)", fit$snp_names)
  dimnames(beta_vals) <- list(varnames, lamNames(fit$lambda))
  
  ## output
  val <- structure(list(beta_vals = beta_vals,
                        lambda = lambda,
                        eta = fit$eta,
                        penalty = fit$penalty,
                        gamma = fit$gamma,
                        alpha = fit$alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = fit$penalty.factor,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        ncol_X = fit$ncol_X,
                        iter = iter,
                        converged = converged),
                   class = "plmm")
  if (fit$returnX) {
    if (utils::object.size(fit$SUX) > 1e8) {
      warning("Due to the large size of SUX (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      # if it fits, it ships 
      returnX <- TRUE
    }
  }
  if(fit$returnX) {
    val$std_X <- fit$std_X # this is the standardized design matrix
    val$y <- fit$y
  } 
  return(val)
  
  
}
    
  