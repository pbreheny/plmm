#' Compute sequence of lambda values
#'
#' This function allows you compute a sequence of lambda values for plmm models.
#' @param X Design matrix which *includes* the intercept column if present. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise. A value of lambda.min = 0 is not supported. 
#' @param nlambda The desired number of lambda values in the sequence to be generated. 
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @export
#' 
#' @examples 
#' RRM <- genRelatednessMat(X = scale(admix$X))
#' fit <- plmm_lasso(X = admix$X, y = admix$y, V = RRM, p1 = 10)
#' (setup_lambda(admix$X, admix$y, alpha = 0.1, nlambda = 10, penalty.factor = fit$fit$penalty.factor)) # use default lambda.min
#' (setup_lambda(admix$X, admix$y, alpha = 0.1, nlambda = 10, penalty.factor = fit$fit$penalty.factor)) # user-specified lambda.min
#' 


setup_lambda <- function(X, y, alpha, lambda.min, nlambda, penalty.factor) {

  # error checking: 
  # make sure alpha is neither missing nor zero
  if(is.na(alpha) | is.null(alpha) | alpha == 0){stop("Must provide a non-zero value for alpha.")}
  # make sure user is not trying to use a lambda.min value of 0
  # if(!missing(lambda.min)){
  #   if(lambda.min == 0){stop("User-specified value for lambda.min cannot be zero")}
  # }
  # 
  
  
  # label dimensions of X 
  n <- nrow(X)
  p <- ncol(X) # including int if intercept = TRUE

  # identify which variables (e.g. SNPs) to penalize 
  ind <- which(penalty.factor != 0)

  # set up a fit from which to derive residuals 
  if (length(ind) != p) { # case 1: not all `p` columns are to be penalized
    fit <- stats::glm(y ~ -1 + X[, -ind, drop = FALSE], family='gaussian')
  } else { # case 2: at least one column is to be penalized 
    fit <- stats::glm(y ~ 1, family='gaussian')
  }

  # determine the 
  zmax <- max(abs((t(X[, ind]) %*% fit$residuals)) / penalty.factor[ind]) / n ### this first part can by xty again
  lambda.max <- zmax/alpha
  # TODO: debugging
  # print(lambda.max)

  # FIXME: working to fix this as of June 23, 2022
  # if (lambda.min == 0) {
  #   lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len = nlambda-1)), 0)
  # } else {
  #   lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len = nlambda))
  # }
  # 
  # if (length(ind)!= p) lambda[1] <- lambda[1] * 1.001
  # lambda
  
  
  # Default is .001 if the number of observations is larger than the number of 
  # covariates and .05 otherwise. A value of lambda.min = 0 is not supported. 
  if(missing(lambda.min)){ # case 1: if the user does not specify a lambda.min value...
    if(n > p){
      
      lambda <- exp(seq(log(lambda.max), log(0.001*lambda.max), len = nlambda))
    } else {
      lambda <- exp(seq(log(lambda.max), log(0.05*lambda.max), len = nlambda))
    }
  } else { # case 2: the user specifies a (nonzero) lambda.min value
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len = nlambda))
  }

}
