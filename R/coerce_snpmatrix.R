#' Coerce SnpMatrix: a function to convert a SnpMatrix object into a matrix of numeric values
#'
#' @param snpmatrix An object of class "SnpMatrix", as defined in the \code{snpStats} package. 
#'
#' @return A matrix of numeric values 
#' @export
#'
#' @examples
#' cad <- process_plink(prefix = "cad_lite", dataDir = plink_example(path="cad_lite.fam", parent=T), coerce=F)
#' X <- coerce_snpmatrix(cad_lite$genotypes)

coerce_snpmatrix <- function(snpmatrix){
  mat <- as(snpmatrix, 'numeric')
  
  # if NAs, throw warning 
  if(sum(is.na(numeric_mat)) > 0) warning("Watch out: there are missing values in the returned matrix.")
  
  return(numeric_mat)
}