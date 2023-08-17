#' std_fbm
#'A function to standardize a matrix stored as an FBM object 
#'
#' @param X 
#' @param ind 
#'
#' @return A standardized matrix of class FBM
#' @keywords internal
std_fbm <- function(X, ind){
ncvreg::std(X[,ind])
}
