#' A summary function for cv.plmm objects
#'
#' @param object A \code{cv.plmm} object
#' @param lambda The regularization parameter value at which inference should be reported. Can choose a numeric value, 'min', or '1se'. Defaults to 'min.'
#' @param ...  Not used 
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
#' @rdname summary.cv.plmm
#' 
#' @export
#'
#' @examples 
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' summary(cv_fit)

summary.cv.plmm <- function(object, lambda = "min", ...){
  
  # determine the number of nonzero coefficients at specified lambda value
  
  if(lambda == "min"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(object$fit, type="nvars", lambda=object$lambda.min)
  }
  
  if(lambda == "1se"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(object$fit, type="nvars", lambda=object$lambda.1se)
  }
  
  if(is.numeric(lambda)){
    if(!(lambda %in% object$fit$lambda)) stop("The user-specified lambda is not one of the lambda values used to fit the model.")
    
    nvars <- predict(object$fit, type="nvars", lambda = lambda)
  }
  
  # TODO: determine what else should go here 
  out <- structure(list(lambda.min = object$lambda.min,
                        lambda.1se = object$lambda.1se,
                        penalty = object$penalty,
                        nvars = nvars,
                        cve = object$cve,
                        min = object$min,
                        fit = object$fit), 
                   class = "summary.cv.plmm")
  if("Bias" %in% names(object)){
    out$bias <- object$Bias
    out$loss <- object$Loss
  }
  
  return(out)
  
}
