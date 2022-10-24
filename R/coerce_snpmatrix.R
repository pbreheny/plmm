#' Coerce SnpMatrix: a function to convert a SnpMatrix object into a matrix of numeric values
#'
#' @param snpmatrix An object of class "SnpMatrix", as defined in the \code{snpStats} package. 
#'
#' @return A matrix of numeric values 
#' @export
#'
#' @examples
#' cad <- process_plink(prefix = "cad", dataDir = plink_example(path="cad.fam", parent=T))
#' X <- coerce_snpmatrix(cad$genotypes)
#' 
coerce_snpmatrix <- function(snpmatrix){
  mat <- as.matrix(snpmatrix)
  numeric_mat <- apply(mat, 2, as.numeric)
  
  # if NAs, throw warning 
  if(sum(is.na(numeric_mat)) > 0) warning("Watch out: there are missing values in the returned matrix.")
  
  return(numeric_mat)
}