#' A function to print the summary of a \code{plmm} model
#'
#' @param x A `summary.plmm` object
#' 
#' @rdname print.summary.plmm
#' @export
#'
#' @examples
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' fit2 <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X), penalty = "SCAD")
#' s1 <- summary(fit, idx = 97)
#' s2 <- summary(fit, lambda = fit$lambda[97])
#' s3 <- summary(fit2, idx = 25)
#' print.summary.plmm(s1)
#' print.summary.plmm(s2)
#' print.summary.plmm(s3)
print.summary.plmm <- function(x){
  
  cat(x$penalty, "-penalized regression model with n=", x$std |> nrow(), ", p=", x$p, sep="")
  cat(" at lambda=", x$lambda_char, "\n", sep="")
  cat("-------------------------------------------------\n")

  # did the model converge? 
  if(x$converged){
      cat("The model converged", "\n")
    } else {
      cat("The model did not converge - max. number of iterations reached", "\n")
    }
 
  cat("-------------------------------------------------\n")
  # nonzero coefficients 
  cat("# of non-zero coefficients: ", x$nvars, "\n")

  # constant features 
  if(length(x$constant_features) > 0){
    if(length(x$constant_features) < 10){
      cat("Constant features (features without variation): ",
          x$constant_features, "\n")
    } else {
      cat("Number of constant features (features without variation): ",
          length(x$constant_features), "\n")
    }
    
  } 
  
  
  cat("-------------------------------------------------\n")
  
  
}
