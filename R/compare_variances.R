#' quick function to do some further simulations 
#'
#' @param nrep # of reps for simulation
#' @param K 
#' @param sig_s 
#' @param sig_eps 
#' @param intercept 
#' @param center_y 
#'
#' @return list with two named vectors: 
#'  * intercept: estimated eta values from the method that defines null model as intercept-only model
#'  * zero_mean: estimated eta values from the model that defines y to have mean 0 on rotated scale. 
#' @export
#'
#' @examples
compare_variances <- function(nrep, K, sig_s, sig_eps,
                              intercept = TRUE, center_y = FALSE){
  hat_eta <- lippert_hat_eta <- rep(NA_integer_, nrep)
  pb <- txtProgressBar(0, nrep, style = 3)
  for(i in 1:nrep){
    # my derivation (has intercept)
    res <- test_eta_estimation(sig_s = sig_s,
                               sig_eps = sig_eps,
                               K = K,
                               intercept = intercept,
                               center_y = center_y)
    hat_eta[i] <- res
    
    # Lippert/Rakitsch way (null model has mean 0)
    lippert_hat_eta[i] <- lippert_test_eta_estimation(sig_s = sig_s,
                                                      sig_eps = sig_eps,
                                                      K = K,
                                                      intercept = intercept,
                                                      center_y = center_y)
    setTxtProgressBar(pb, i)
  }
  
  return(list(
    intercept = hat_eta,
    zero_mean = lippert_hat_eta
  ))
  
}
