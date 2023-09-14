#' Predict method for a list used in cross-validation (e.g. the \code{fit.i} within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`
#' @param newX A design matrix used for computing predicted values (i.e, the test data).
#' @param type A character argument indicating what type of prediction should be returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#' @param prep Optional argument. Result of the call to `plmm_prep` which corresponds to the `fit` argument. Required if \code{type == 'blup'}. 
#' @param V11 Variance-covariance matrix of the training data. Extracted from `estimated_V` that is generated using all observations. Required if \code{type == 'blup'}. 
#' @param V21 Covariance matrix between the training and the testing data. Extracted from `estimated_V` that is generated using all observations. Required if \code{type == 'blup'}. 
#' @param ... Additional optional arguments
#' 
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows: 
#'  * 'lp' (default): uses the linear predictor (i.e., product of new data and estimated coefficients) to predict new values of the outcome. 
#'    Note that this approach does not incorporate the correlation structure of the data. 
#'  
#'  * 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'lp' a value that represents the estimated random effect. 
#'  This addition is a way of incorporating the estimated correlation structure of data into our prediction of the outcome. 
#'  
#' @keywords internal
#'

predict.list <- function(fit, newX, type=c("lp", "blup"),
                         lambda, idx=1:length(fit$lambda), prep = NULL, V11 = NULL, V21 = NULL, ...) {
  type <- match.arg(type)
  beta_vals <- coef.list(fit, lambda=lambda, which=idx, drop=FALSE) # includes intercept 
  
  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)
  
  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  Xbeta <- cbind(1, newX) %*% beta_vals
  
  if (type=="lp") return(drop(Xbeta))
  
  if (type == "blup"){
    # warning("The BLUP option is under development. Rely on these estimates at your own risk.")
    
    if(is.null(prep)) stop("The 'prep' argument is required for BLUP calculation.")
      
    # covariance comes from selected rows and columns from estimated_V that is generated in the overall fit (V11, V21)
      
    # test1 <- V21 %*% chol2inv(chol(V11)) # true 
    # TODO: to find the inverse of V11 using svd results of K, i.e., the inverse of a submatrix, might need to use Woodbury's formula 
    
    std_y <- ncvreg::std(fit$y) # get y on same scale as X!
    ranef <- V21 %*% chol2inv(chol(V11)) %*% (std_y - cbind(1, prep$std_X) %*% beta_vals)

      
    blup <- Xbeta + ranef
    
    return(blup)
  }
  
}



