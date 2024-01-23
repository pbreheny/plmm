#' Read in processed data
#' This function is intended to be called after `process_plink()` has been called once. 
#' 
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path. 
#' @param row_id Character string indicating which IDs to use for the rownames of the genotype matrix. Can choose "fid" or "iid", corresponding to the first or second columns in the PLINK .fam file. Defaults to NULL. 
#' @param returnX Logical: Should the design matrix be returned as a numeric matrix that will be stored in memory. By default, this will be FALSE if the object sizes exceeds 100 Mb.
#' @param trace Logical: Should trace messages be shown? Default is TRUE. 
#' 
#' @return A list with these components: 
#'  * std_X, the column-standardized design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details. 
#'  * fam, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * map, a data frame containing the feature information (like a .bim file in PLINK)
#'  * ns: A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0. 
#'  * center: A vector of values for centering each column in X
#'  * scale: A vector of values for scaling each column in X 
#' 
#' @importFrom data.table setorderv
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' pen <- get_data(path = "../temp_files/penncath_lite", trace = TRUE)
#' }
#' 
#' @details
#' The .rds object should have an 'std_X' element - this is what will be used as the design matrix for analysis. This design matrix should *not* include an intercept column (this will be added later in `plmm_fit`()).
#' 
#' In the returned list, the `fam` data will be sorted by family and by individual, as in `dplyr::arrange(family.ID, sample.ID)`.
#' The rows of `X` will be sorted to align in the same order as in `fam`, where rownames of `X` will be sample ID. 
#' 
#' 
get_data <- function(path, row_id = NULL, returnX, trace = TRUE){
  
  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods 
  
  if(file.exists(bk)){
    obj <- bigsnpr::snp_attach(rdsfile = rds)
  } else {
    obj <- readRDS(rds)
  }
  
  # return data in a tractable format 
  if (missing(returnX)) {
    if (utils::object.size(obj$std_X) > 1e8) {
      warning("\nDue to the large size of X (>100 Mb), X has been returned as a file-backed matrix (FBM).
              \nTo turn this message off, explicitly specify fbm=TRUE or fbm=FALSE).")
      returnX <- TRUE
    } else {
      # if it fits, it ships 
      returnX <- FALSE
    }
  }
  
  if(returnX){
    # get std_X as a matrix 
    std_X <- obj$std_X[,]
    if(!is.null(row_id)){
      if(row_id == "iid"){row_names <- obj$fam$sample.ID}
      if(row_id == "fid"){row_names <- obj$fam$family.ID}
    } else {
      row_names <- 1:length(obj$fam$sample.ID)
    }
    
    dimnames(std_X) <- list(row_names,
                        obj$map$marker.ID[obj$ns])
    # TODO fix this error: Error in dimnames(X) <- list(row_names, o
    
    
    if(!(all.equal(obj$fam$sample.ID, as.numeric(rownames(std_X))))){
      stop("\nThere is an issue with the alignment between the rownames of the genotype data and the sample IDs.
           \nWere there individuals represented in the .bed file who are not in the .fam file, or vice versa?
           \nPlease ensure that your PLINK files represent all the same individuals before analyzing data with PLMM.")
    }
    
    
    cat("\nReminder: the X that is returned here is column-standardized, with constant features removed.
        \nA copy of the original data is available via the 'genotypes' matrix in the .rds object")
    return(list(n = obj$n,
                p = obj$p,
                std_X = std_X,
                std_X_center = obj$std_X_center,
                std_X_scale = obj$std_X_scale,
                fam = obj$fam,
                map = obj$map,
                ns = obj$ns))
  } else {
    cat("Note: The design matrix is being returned as a file-backed matrix (FBM) -- see bigstatsr::FBM() for details.")
    
    cat("\nReminder: the X that is returned here is column-standardized.
        \nA copy of the original data is available via the 'genotypes' matrix in the .rds object")
    
    return(list(n = obj$n,
                p = obj$p,
                std_X = obj$std_X,
                std_X_center = obj$std_X_center,
                std_X_scale = obj$std_X_scale,
                fam = obj$fam,
                map = obj$map,
                ns = obj$ns))
  }
  
}