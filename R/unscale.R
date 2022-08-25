#' Unscale coefficient values
#'
#' This function allows you to unscale coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Unscaled, rotated design matrix *with* an intercept column if present. Necessary for properly defining the dimensions of beta in cases where singular columns are present.
#' @param scaled_X Scaled, rotated design matrix. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export
#'
#' @examples 
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' unscaled_betas <- unscale(b = fit$beta, X = admix$X, scaled_X = fit$SUX, intercept = TRUE)

unscale <- function(b, X, scaled_X, intercept = TRUE) {
  # identify which columns of the scaled, rotated design matrix are nonsingular 
  ns <- attr(scaled_X, 'nonsingular')
  # extract the scaling values from the non-singular columns
  scale <- attr(scaled_X, 'scale')[ns]
  
  # calculate unscaled beta values 
  if (intercept){ # case 1: intercept 
    # create a matrix for beta values to be "filled in"
    unscaled_beta <- matrix(0, nrow = ncol(X) + 1, ncol = ncol(b))
    
    unscaled_beta[1,] <- b[1, , drop = FALSE] # don't change the intercept values 
    b <- b[-1, , drop=FALSE] # now, b has only non-intercept values 
    bb <- b / scale # unscale the non-intercept values 
    unscaled_beta[1 + ns,] <- bb # fill in the rest of the beta values 
  } else { # case 2: no intercept 
    # create a matrix for beta values to be "filled in"
    unscaled_beta <- matrix(0, nrow = ncol(X), ncol = ncol(b))
    bb <- b / scale # unscale the coefficient values 
    unscaled_beta[ns, ] <- bb # fill in beta values for nonsingular coefficients 
  }
  return(unscaled_beta)
}
