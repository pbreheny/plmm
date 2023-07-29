#' Read in processed data
#' This function is intended to be called after `process_plink()` has been called once. 
#' 
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path. 
#' @param row_id Character string indicating which IDs to use for the rownames of the genotype matrix. Can choose "fid" or "iid", corresponding to the first or second columns in the PLINK .fam file. Defaults to NULL. 
#' @param returnX Logical: Should the design matrix be returned as a numeric matrix that will be stored in memory. By default, this will be FALSE if the object sizes exceeds 100 Mb.
#' @param trace Logical: Should trace messages be shown? Default is TRUE. 
#' 
#' @return A list with four components: 
#'  * X, the design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details. 
#'  * fam, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * map, a data frame containing the feature information (like a .bim file in PLINK)
#'  * constants_idx A vector indicating the which columns of X contain constant features (features with no variance). 
#' @export
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
    if (utils::object.size(obj$genotypes) > 1e8) {
      warning("\nDue to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      # if it fits, it ships 
      returnX <- TRUE
    }
  }
  
  if(returnX){
    X <- obj$genotypes[,]
    if(!is.null(row_id)){
      if(row_id == "iid"){row_names <- obj$fam$sample.ID}
      if(row_id == "fid"){row_names <- obj$fam$family.ID}
    } else {
      row_names <- 1:length(obj$fam$sample.ID)
    }
    
    dimnames(X) <- list(row_names,
                        obj$map$marker.ID)
    
    return(list(X = X,
                fam = obj$fam,
                map = obj$map,
                constants_idx = which(obj$constants_idx)
                ))
  } else {
    return(list(X = obj$genotypes,
                fam = obj$fam,
                map = obj$map,
                constants_idx = which(obj$constants_idx)
    ))
  }
  
}