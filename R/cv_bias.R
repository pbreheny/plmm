#' Bias in the i-th fold in cross-validation due to correlation among observations w_cv_GLS on page 722. This does not depend on lambda  
#'
#' @param Xtrain Training set design matrix 
#' @param Xtest Testing set design matrix 
#' @param ytrain Training set response 
#' @param ytest Testing set response 
#' @param cov_ytrain Covariance matrix of ytrain 
#' @param cov_ytrain_ytest Covariance matrix between ytrain and ytest 
#' 
#' 
#' @export
#' 
#' @examples 
#' 





cv_bias <- function(Xtrain, Xtest, ytrain, ytest, cov_ytrain, cov_ytrain_ytest) {
  cov_ytrain_inv <- chol2inv(chol(cov_ytrain)) 
  Xtrain_cov_ytrain_inv <- crossprod(Xtrain, cov_ytrain_inv) 
  first_half <- crossprod(Xtest, chol2inv(chol(Xtrain_cov_ytrain_inv %*% Xtrain))) 
  second_half <- Xtrain_cov_ytrain_inv %*% cov_ytrain_ytest
  return(first_half * second_half)
}


# Notes: 
# S is D^2 divided by p or not? In Anna's thesis page 18, it's not divided by p 
# V^{-1/2} = UW, which does not depend on lambda. 
## Let plmm_fit return estimated_cov_y , to be used in cv_bias line 34 
# line 129 in cv_plmm, is that w_cv_GLS from the paper? 


