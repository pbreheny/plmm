#' Predict method for a list used in cross-validation (within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`. 
#' @param oldX The standardized design matrix of training data, *pre-rotation*. 
#' @param newX A design matrix used for computing predicted values (i.e, the test data).
#' @param type A character argument indicating what type of prediction should be returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
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
                         oldX,
                         newX,
                         type=c("lp", "blup"),
                         idx=1:length(fit$lambda),
                         V11 = NULL,
                         V21 = NULL, ...) {
  
  type <- match.arg(type)
  # get beta values (for nonsingular features) from fit
  beta_vals <- fit$untransformed_b1[,idx,drop = FALSE]
  
  # format dim. names
  if(is.null(dim(beta_vals))) {
    # case 1: beta_vals is a vector 
    names(beta_vals) <- lamNames(fit$lambda)
  } else {
    # case 2: beta_vals is a matrix
    colnames(beta_vals) <- lamNames(fit$lambda)
  }
  # calculate the estimated mean values for test data 
  a <- beta_vals[1,]
  b <- beta_vals[-1,,drop=FALSE]
 
  Xb <- sweep(newX %*% b, 2, a, "+")
  
  # for linear predictor, return mean values 
  if (type=="lp") return(drop(Xb))
  
  # for blup, will incorporate the estimated variance 
  if (type == "blup"){
    # covariance comes from selected rows and columns from estimated_V that is generated in the overall fit (V11, V21)
    # test1 <- V21 %*% chol2inv(chol(V11)) # true 
    # TODO: to find the inverse of V11 using svd results of K, i.e., the inverse of a submatrix, might need to use Woodbury's formula 
    Xb_train <- sweep(oldX %*% b, 2, a, "+")
    resid_train <- (fit$y - Xb_train)
    ranef <- V21 %*% (chol2inv(chol(V11)) %*% resid_train)
    blup <- Xb + ranef
    
    return(blup)
  }
  
}



