#' A function to help with accessing example PLINK files 
#'
#' @param path Optional argument specifying a path (filename) for an external data file in \code{inst/extdata} 
#' @param parent If \code{path=TRUE} and the user wants the name of the parent directory where that file is located, set \code{parent=TRUE}. Defaults to FALSE. 
#'
#' @return If \code{path=NULL}, a character vector of file names is returned. If path is given, then a character string 
#' with the full file path
#' 
#' @export
#'
#' @examples
#' # name of all external PLINK data files 
#' (plink_example())
#' 
#' # name of file path for a specific PLINK file 
#' plink_example(path="cad.fam") # what shows up here should be unique to the user's machine 
#' 
#' # name of file path for parent directory of a specific file 
#' plink_example(parent=T, path="cad.fam")

plink_example <- function(path = NULL, parent=FALSE) {
  if (is.null(path)) {
    dir(base::system.file("inst/extdata", package = "penalizedLMM"))
  } else {
    if(parent){
      base::system.file("inst/extdata", package = "penalizedLMM", mustWork = T)
    } else {
      base::system.file("inst/extdata", path, package = "penalizedLMM", mustWork = TRUE)
    }
  }
}
