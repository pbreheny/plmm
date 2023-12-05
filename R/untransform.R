#' Untransform coefficient values back to the original scale 
#' 
#' This function unwinds both steps of the standardization process to obtain
#' coefficient values on their original scale. It is called by plmm_format().
#' 
#' @param res_b The values returned in the 'beta' argument of the ncvfit() object
#' @param ns The indices of the non-singular columns of the ORIGINAL design matrix
#' @param p The number of columns in the original design matrix (without the intercept)
#' @param std_X_details A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements 'scale', 'center', and 'nonsingular'
#' @param rot_X The rotated design matrix (has the intercept)
#' @param stdrot_X The standardized design matrix AFTER rotation; this should have attribute 'scale'
#' @param partial Logical: Should beta values be untransformed from just the post-rotation standardization? 
#' This option exists because I need partial= TRUE inside the internal function `predict.list()`. 
#' Users should leave this option false unless they know exactly what they're doing. 
#' @keywords internal


untransform <- function(res_b, ns, p, std_X_details, rot_X, stdrot_X, partial = FALSE){
  ## Step 1: reverse the POST-ROTATION standardization ##  

  # create placeholder vector
  untransformed_b1 <- res_b
  
  # un-scale the non-intercept values & fill in the placeholder
  untransformed_b1[-1,] <- sweep(x = res_b[-1, , drop=FALSE],
                                 MARGIN = 1, # beta values are on rows 
                                 STATS = attr(stdrot_X, 'scale'),
                                 FUN = "/")
  
  if (partial){
    return(untransformed_b1)
  } else {
    ## Step 2: reverse the PRE-ROTATION standardization ## 
    
    # partition the values from Step 1 into intercept and non-intercept parts
    a <- untransformed_b1[1, , drop = FALSE] # this is the intercept 
    b <- untransformed_b1[-1, , drop=FALSE]
    
    # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
    # this will create columns of zeros for betas corresponding to singular columns
    untransformed_beta <- matrix(0,
                                 nrow = (p + 1), # + 1 is for the intercept
                                 ncol = ncol(untransformed_b1))

    # unscale the beta values for non-singular, non-intercept columns
    untransformed_b2 <- sweep(x = b,
                              MARGIN = 1,
                              STATS = std_X_details$scale[ns],
                              FUN = "/")
    
    # fill in the un-transformed values
    untransformed_beta[ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
    untransformed_beta[1,] <- a - crossprod(std_X_details$center[ns], untransformed_b2)

    
    # Final step: return un-transformed beta values 
    return(untransformed_beta)
  }
  
}


