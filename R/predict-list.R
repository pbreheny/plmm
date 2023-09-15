#' Predict method for a list used in cross-validation (within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`. 
#' @param std_X The standardized design matrix of training data, *pre-rotation*. 
#' @param newX A design matrix used for computing predicted values (i.e, the test data).
#' @param type A character argument indicating what type of prediction should be returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
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

predict.list <- function(fit,
                         std_X,
                         newX,
                         type=c("lp", "blup"),
                         idx=1:length(fit$lambda),
                         prep = NULL,
                         V11 = NULL,
                         V21 = NULL, ...) {
  
  type <- match.arg(type)
  # get beta values from fit -- these coefficients are on the transformed scale
  raw_beta_vals <- fit$b[, idx, drop = FALSE]
  
  # reverse the transformations of the beta values 
  beta_vals <- untransform(res_b = raw_beta_vals,
                           ns = fit$ns,
                           ncol_X = fit$ncol_X,
                           std_X = std_X,
                           SUX = fit$SUX,
                           std_SUX = fit$std_SUX,
                           partial = TRUE)
  
  # format dim. names
  if(is.null(dim(beta_vals))) {
    # case 1: beta_vals is a vector 
    names(beta_vals) <- lamNames(fit$lambda)
  } else {
    # case 2: beta_vals is a matrix
    colnames(beta_vals) <- lamNames(fit$lambda)
  }

  # calculate the estimated mean values for test data 
  Xbeta <- cbind(1, newX) %*% beta_vals
  
  # for linear predictor, return mean values 
  if (type=="lp") return(drop(Xbeta))
  
  # for blup, will incorporate the estimated variance 
  if (type == "blup"){
    # warning("The BLUP option is under development. Rely on these estimates at your own risk.")
    
    if(is.null(prep)) stop("The 'prep' argument is required for BLUP calculation.")
      
    # covariance comes from selected rows and columns from estimated_V that is generated in the overall fit (V11, V21)
      
    # test1 <- V21 %*% chol2inv(chol(V11)) # true 
    # TODO: to find the inverse of V11 using svd results of K, i.e., the inverse of a submatrix, might need to use Woodbury's formula 

    ranef <- V21 %*% chol2inv(chol(V11)) %*% (fit$y - cbind(1, prep$std_X) %*% beta_vals)

    blup <- Xbeta + ranef
    
    return(blup)
  }
  
}



