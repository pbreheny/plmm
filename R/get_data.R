#' Read in processed data
#' This function is intended to be called after `process_plink()` has been called once. 
#' 
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path. 
#' @param returnX Logical: Should the design matrix be returned as a numeric matrix that will be stored in memory. By default, this will be FALSE if the object sizes exceeds 100 Mb.
#' @param trace Logical: Should trace messages be shown? Default is TRUE. 
#' 
#' @returns A list with four components: 
#'  * `X`, the design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details. 
#'  * `fam`, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * `map`, a data frame containing the feature information (like a .bim file in PLINK)
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' pen <- get_data(path = "../temp_files/penncath_lite", trace = TRUE)
#' }
#' 
#' @details
#' In the returned list, the `fam` data will be sorted by family and by individual, as in `dplyr::arrange(family.ID, sample.ID)`.
#' The rows of `X` will be sorted to align in the same order as in `fam`, where rownames of `X` will be sample ID. 
#' 
#' 
get_data <- function(path, returnX, trace = TRUE){
  
  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods 
  
  if(file.exists(bk)){
    obj <- bigsnpr::snp_attach(rdsfile = rds)
  } else {
    obj <- readRDS(rds)
  }
  
  # return data in a tractable format 
  if (missing(returnX)) {
    if (utils::object.size(obj$genotypes) > 1e8) {
      warning("\nDue to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      # if it fits, it ships 
      returnX <- TRUE
    }
  }
  
  # TODO: should I order fam data by FID and IID? 

  if(returnX){
    X <- obj$genotypes[,]
    dimnames(X) <- list(obj$fam$sample.ID, obj$map$marker.ID)
    
    # TODO: think about ordering the rows of X to match fam file. 
    # Is this ordering something that is wise to do? 
    
    # if(!(all.equal(obj$fam$sample.ID, as.numeric(rownames(X))))){
    #   stop("\nThere is an issue with the alignment between the rownames of the genotype data and the sample IDs. 
    #        \nWere there individuals represented in the .bed file who are not in the .fam file, or vice versa?
    #        \nPlease ensure that your PLINK files represent all the same individuals before analyzing data with PLMM.")
    # }
    return(list(X = X,
                fam = obj$fam,
                map = obj$map))
  } else {
    cat("Note: X is being returned as a file-backed matrix (FBM) -- see bigstatsr::FBM() for details.
        \n At this time, plmm() cannot analyze design matrix X in this FBM format. Allowing such an option will
        \n require writing the 'meat and potatoes' of plmm() in C++, which is a work in progress. For now, 
        \n functions from package bigsnpr may be used for analyzing FBM data.")
    
    
    # warning("\nSorting the returned X matrix is not yet implemented for file-backed data.
    #         \nFor analyses, remember to sort to sort the rows of the genotype data 
    #         in the same order as the rows of the family data, so that matrices X and Z have same pattern!")
    
    return(list(X = obj$genotypes,
                fam = obj$fam,
                map = obj$map))
  }
  
}