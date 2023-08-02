#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param convex convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. 
#' @param K Similarity matrix used to rotate the data. This will be passed from plmm() as (1) NULL, in which case relatedness_mat(std(X)) is used, (2) a matrix or (3) a list with components d and u (eigenvalues and eigenvectors, respectively).
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
#' @keywords internal
#'
#' @examples
#' 
#' \dontrun{
#' # NB: this is an internal function (i.e. it's not exported)
#' # these examples won't work unless you have penalizedLMM::: 
#' prep1 <- plmm_prep(X = admix$X, y = admix$y, trace = TRUE)
#' fit1 <- plmm_fit(prep = prep1)
#' res1 <- plmm_format(fit1)
#' 
#' }

plmm_format <- function(fit,
                        convex = TRUE,
                        dfmax = fit$ncol_X + 1, 
                        X,
                        K){
  
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
  # browser()
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
  
  if (is.null(K)) {
    Vhat <- fit$eta * (1/fit$ncol_X) * tcrossprod(X) + (1-fit$eta) * diag(fit$nrow_X) 
  } else if (is.list(K)){
    K <- fit$U %*% tcrossprod(diag(fit$S), fit$U)
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$nrow_X)
  } else if (is.matrix(K)){
    Vhat <- fit$eta * K + (1-fit$eta) * diag(fit$nrow_X)
  }
  ## output
  val <- structure(list(beta_vals = beta_vals,
                        lambda = lambda,
                        eta = fit$eta,
                        S = fit$S,
                        U = fit$U,
                        penalty = fit$penalty,
                        gamma = fit$gamma,
                        alpha = fit$alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = fit$penalty.factor,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        ncol_X = fit$ncol_X,
                        nrow_X = fit$nrow_X, 
                        Vhat = Vhat, 
                        iter = iter,
                        converged = converged),
                   class = "plmm")
  if (fit$returnX) {
    if (utils::object.size(fit$SUX) > 1e8) {
      warning("Be aware that SUX is large (>100 Mb)")
    }
    val$std_X <- fit$std_X # this is the standardized design matrix
    val$y <- fit$y
  }
  
  if("S" %in% names(fit) & "U" %in% names(fit)){
    val$S <- fit$S
    val$U <- fit$U
  }
  
  return(val)
  
  
}
    
  