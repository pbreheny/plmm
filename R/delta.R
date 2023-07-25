#' a helper function to calculate the 'distance' between 2 matrices 
#' @param A The matrix to be compared with B
#' @param B The second matrix (I think of this as the approximation to A)
#' @param type String indicating the type of norm to be used. See Matrix::norm() for details
#' 
#' @keywords internal
#' 
delta <- function(A, B, type = "F"){
  Matrix::norm(A, type = type) - Matrix::norm(B, type = type)
}


