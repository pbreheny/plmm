#' Predict method for a list used in cross-validation (e.g. the \code{fit.i} within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`
#' @param newX A design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what type of prediction should be returned. Options are "response," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#' @param prep Optional argument. Result of the call to `plmm_prep` which corresponds to the `fit` argument. Required if \code{type == 'blup'} and object is too large to be returned in `fit` object.
#' @param covariance Optional argument. $q times n$ covariance matrix between new and old observations. Required if \code{type == 'blup'}.
#' @param ... Additional optional arguments
#' @export
#'

predict.list <- function(fit, newX, type=c("response", "coefficients", "vars", "nvars", "blup"),
                         lambda, idx=1:length(fit$lambda), prep = NULL, ...) {
  type <- match.arg(type)
  beta_vals <- coef.list(fit, lambda=lambda, which=idx, drop=FALSE) # includes intercept 
  
  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)
  
  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  Xbeta <- cbind(1, newX) %*% beta_vals
  
  if (type=="response") return(drop(Xbeta))
  
  if (type == "blup"){
    warning("The BLUP option is under development. Rely on these estimates at your own risk.")
    
      if(is.null(prep)) stop("The 'prep' argument is required for BLUP calculation.")
      
      # calculate covariance between new and old observations 
     # TODO: ask Patrick if using standardized design matrices in the 
    # covariance calculation below is appropriate.... 
    
      covariance <- cov(t(ncvreg::std(newX)), t(prep$std_X))
      
      ranef <- covariance %*% prep$U %*% diag((1 + fit$eta * (prep$S - 1))^(-1)) %*% t(prep$U) %*% (fit$y - cbind(1, prep$std_X) %*% beta_vals)
      # print(eta) 
      
      blup <- Xbeta + ranef
      
    
    return(blup)
  }
  
}



