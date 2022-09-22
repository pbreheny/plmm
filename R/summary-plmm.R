#' A summary method for the plmm class 
#'
#' @param object An object of class \code{plmm}
#' @param quiet A logical indicating whether the console output is desired. Defaults to FALSE.
#' @return A list containing: 
#' @export
#'
#' @examples
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' example <- summary.plmm(fit)
#' #TODO: Sept. 14, 2022 - need a more sophisticated way to select the 'best' lambda value 

summary.plmm <- function(object, quiet = FALSE){
  
  # tell the chosen lambda value 
  best_lam_idx <- which.min(object$loss)
  best_lam <- object$lambda[best_lam_idx]
  if(!quiet){
    cat("The lambda value that minimizes loss: ", best_lam, "\n")
  }
  
  
  # did the model converge? 
  if(!quiet){
    if(object$converged[best_lam_idx]){
      cat("The model converged", "\n")
    } else {
      cat("The model did not converge - max. number of iterations reached", "\n")
    }
  }
  
  # tell the number of betas estimated 
  if(!quiet){
    cat("# of coefficients estimated: ", dim(object$beta_vals)[1], "\n") 
  }
  
  
  # tell the number of non-zero betas 
  select_beta_vals <- object$beta_vals[,best_lam_idx]
  nz_beta_vals <- select_beta_vals[select_beta_vals > .Machine$double.eps]

  if(!quiet){
    cat("# of non-zero coefficients: ", length(nz_beta_vals), "\n")
  }
  
  
  # list the predictors corresponding to the non-zero betas 
  nz_beta_idx <- which(select_beta_vals %in% nz_beta_vals)
  chosen_predictors <- rownames(object$beta_vals)[nz_beta_idx]
  
  if(!quiet){
    cat("The predictors included in the model are: ", chosen_predictors, "\n")
  }
  
  
  # if there were any monomorphic SNPs, say so: 
  monomorphic_snps <- rownames(object$beta_vals)[rowSums(object$beta_vals) == 0]
  
  if(!quiet){
    if(length(monomorphic_snps) == 0){
      cat("There were no monomorphic SNPs", "\n")
    } else {
      cat("The monomorphic SNPs were: ", monomorphic_snps, "\n")
    }
  }
  
  
  return(list(
    best_lam_idx = best_lam_idx,
    best_lam = best_lam,
    select_beta_vals = select_beta_vals,
    chosen_predictors = chosen_predictors,
    monomorphic_snps = monomorphic_snps
  ))
  
}

