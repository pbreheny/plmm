#' A summary function for cv.plmm objects
#'
#' @param obj A \code{cv.plmm} object
#' @param lambda The regularization parameter value at which inference should be reported. Can choose a numeric value, 'min', or '1se'. Defaults to 'min.'
#' @return The return value is an object with S3 class `summary.cv.plmm`. The class has its own print method and contains the following list elements: 
#' * `lambda.min`: The lambda value at the minimum cross validation error 
#' * `lambda.1se`: The maximum lambda value within 1 standard error of the minimum cross validation error 
#' * `nvars`: The number of non-zero coefficients at the selected lambda value 
#' * `cve`: The cross validation error at all folds
#' * `min`: The minimum cross validation error 
#' * `fit`: The \code{plmm} fit used in the cross validation
#' * `bias`: The mean bias of the cross validation
#' * `loss`: The loss (at each fold?) TODO: double-check this  
#' * `penalty`: The penalty applied to the fitted model
#' 
#' 
#' @export
#'
#' @examples 
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' s <- summary.cv.plmm(cv_fit)
#' print.summary.cv.plmm(s)
summary.cv.plmm <- function(obj, lambda = "min"){
  
  # determine the number of nonzero coefficients at specified lambda value
  
  if(lambda == "min"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(obj$fit, type="nvars", lambda=obj$lambda.min)
  }
  
  if(lambda == "1se"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(obj$fit, type="nvars", lambda=obj$lambda.1se)
  }
  
  if(is.numeric(lambda)){
    if(!(lambda %in% obj$fit$lambda)) stop("The user-specified lambda is not one of the lambda values used to fit the model.")
    
    nvars <- predict(obj$fit, type="nvars", lambda = lambda)
  }
  
  # TODO: determine what else should go here 
  out <- structure(list(lambda.min = obj$lambda.min,
                        lambda.1se = obj$lambda.1se,
                        penalty = obj$penalty,
                        nvars = nvars,
                        cve = obj$cve,
                        min = obj$min,
                        fit = obj$fit,
                        bias = obj$Bias,
                        loss = obj$Loss), 
                   class = "summary.cv.plmm")
  return(out)
  
}
