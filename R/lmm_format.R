#' LMM format: a function to format the output of a model constructed with \code{lmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{lmm_fit}
#' @param convex convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. 
#' @param K Similarity matrix used to rotate the data. This will be passed from plmm() as (1) NULL, in which case relatedness_mat(std(X)) is used, (2) a matrix or (3) a list with components d and u (eigenvalues and eigenvectors, respectively).
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
                        convex = TRUE,
                        dfmax = fit$ncol_X + 1, 
                        X,
                        K){
  
  if (fit$warn & sum(iter) == fit$max.iter) warning("Maximum number of iterations reached")
  # TODO: determine if I should just take out the lines below 
  # convex.min <- if (convex) convexMin(b = fit$b,
  #                                     X = fit$std_SUX,
  #                                     family = 'gaussian') else NULL
  
  # reverse the transformations of the beta values 
  beta_vals <- untransform(res_b = fit$b,
                           ns = fit$ns,
                           ncol_X = fit$ncol_X,
                           std_X = fit$std_X,
                           SUX = fit$SUX,
                           std_SUX = fit$std_SUX)
  
  if(fit$trace){cat("\nBeta values are estimated -- almost done!\n")}
  
  # give the matrix of beta_values readable names 
  # SNPs (or covariates) on the rows
  varnames <- c("(Intercept)", fit$snp_names)
  dimnames(beta_vals) <- list(varnames, 'Coeff.')
  
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
                        eta = fit$eta,
                        SUX = fit$SUX,
                        SUy = fit$SUy,
                        convex.min = convex.min, # TODO: may need to take out convex.min
                        loss = loss,
                        penalty.factor = fit$penalty.factor,
                        ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE! 
                        ncol_X = fit$ncol_X,
                        nrow_X = fit$nrow_X, 
                        # estimated_V = fit$estimated_V, # this is using the standardized X scale, used in cv-plmm 
                        Vhat = Vhat, 
                        iter = iter,
                        converged = converged),
                   class = "lmm")
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
  
  if("S" %in% names(fit) & "U" %in% names(fit)){
    val$S <- fit$S
    val$U <- fit$U
  }
  
  return(val)
  
  
}

