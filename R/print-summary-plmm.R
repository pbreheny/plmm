#' A function to print the summary of a \code{plmm} model
#'
#' @param x A `summary.plmm` object
#' @param ... Not used
#' @rdname print.summary.plmm
#' @export
#'
#' @examples
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' fit2 <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X), penalty = "SCAD")
#' s1 <- summary(fit, idx = 97)
#' s2 <- summary(fit, lambda = fit$lambda[97])
#' s3 <- summary(fit2, idx = 25)
#' print(s1)
#' print(s2)
#' print(s3)
print.summary.plmm <- function(x, ...){
  
  cat(x$penalty, "-penalized regression model with n=", x$n, ", p=", x$p, sep="")
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
  cat("-------------------------------------------------\n")
  
  
}
