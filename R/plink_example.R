#' A function to help with accessing example PLINK files 
#'
#' @param path Argument (string) specifying a path (filename) for an external data file in \code{extdata/} 
#' @param parent If \code{path=TRUE} and the user wants the name of the parent directory where that file is located, set \code{parent=TRUE}. Defaults to FALSE. 
#'
#' @return If \code{path=NULL}, a character vector of file names is returned. If path is given, then a character string 
#' with the full file path
#' 
#' @export
#'

plink_example <- function(path, parent=FALSE) {
  if(parent){
    # if parent option selected, will return path to folder 
    system.file("extdata", package = "penalizedLMM", mustWork = TRUE)
  } else {
    # if parent option not selected, will return path to folder WITH file name
    system.file("extdata", path, package = "penalizedLMM", mustWork = TRUE)
  }
}
