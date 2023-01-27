#' Untransform coefficient values back to the original scale 
#' 
#' This function unwinds both steps of the standardization process to obtain
#' coefficient values on their original scale. It is called by plmm().
#' 
#' @param res_b The values returned in the 'beta' argument of the ncvfit() object
#' @param ns The indices of the non-singular columns of the ORIGINAL design matrix
#' @param ncol_X The number of columns in the original design matrix (without the intercept)
#' @param std_X The standardized design matrix BEFORE rotation; this should have attributes 'scale', 'center', and 'nonsingular'
#' @param SUX The rotated design matrix (has the intercept)
#' @param std_SUX The standardized design matrix AFTER rotation; this should have attribute 'scale'
#' 
#' @keywords internal


untransform <- function(res_b, ns, ncol_X, std_X, SUX, std_SUX){
  
## Step 1: reverse the POST-ROTATION standardization ##  
  
  # create placeholder vector
  untransformed_b1 <- res_b
  
  # un-scale the non-intercept values & fill in the placeholder
  untransformed_b1[-1,] <- res_b[-1, , drop=FALSE]/attr(std_SUX, 'scale')
  
## Step 2: reverse the PRE-ROTATION standardization ## 
  
  # partition the values from Step 1 into intercept and non-intercept parts
  a <- untransformed_b1[1, , drop = FALSE] # this is the intercept 
  b <- untransformed_b1[-1, , drop=FALSE]
  
  # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
  # this will create columns of zeros for betas corresponding to singular columns
  untransformed_beta <- matrix(0,
                               nrow = (ncol_X + 1), # + 1 is for the intercept
                               ncol = ncol(untransformed_b1))
  
  # unscale the beta values for non-singular, non-intercept columns
  untransformed_b2 <- b/attr(std_X, 'scale')[ns]
  
  # fill in the un-transformed values
  untransformed_beta[ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
  untransformed_beta[1,] <- a - crossprod(attr(std_X, 'center')[ns], untransformed_b2)
  
  # Final step: return un-transformed beta values 
  return(untransformed_beta)
}


