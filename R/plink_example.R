#' A function to help with accessing example PLINK files 
#'
#' @param path Optional argument specifying a path (filename) for an external data file in \code{extdata/} 
#' @param parent If \code{path=TRUE} and the user wants the name of the parent directory where that file is located, set \code{parent=TRUE}. Defaults to FALSE. 
#'
#' @return If \code{path=NULL}, a character vector of file names is returned. If path is given, then a character string 
#' with the full file path
#' 
#' @export
#'

plink_example <- function(path = NULL, parent=FALSE) {
  if (is.null(path)) {
    dir(base::system.file("extdata", package = "penalizedLMM"))
  } else {
    if(parent){
      base::system.file("extdata", package = "penalizedLMM", mustWork = TRUE)
    } else {
      base::system.file("extdata", path, package = "penalizedLMM", mustWork = TRUE)
    }
  }
}
