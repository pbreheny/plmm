#' Read in processed data
#' This function is intended to be called after `process_plink()` has been called once. 
#' 
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path. 
#' @param row_id Character string indicating which IDs to use for the rownames of the genotype matrix. Can choose "fid" or "iid", corresponding to the first or second columns in the PLINK .fam file. Defaults to NULL. 
#' @param fbm Logical: Should the design matrix be returned as an object of type Filebacked Big Matrix (FBM) as opposed to a numeric matrix that will be stored in memory. By default, this will be TRUE if the object sizes exceeds 100 Mb.
#' @param trace Logical: Should trace messages be shown? Default is TRUE. 
#' 
#' @return A list with four components: 
#'  * std_X, the column-standardized design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details. 
#'  * fam, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * map, a data frame containing the feature information (like a .bim file in PLINK)
#'  * ns: A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0. 
#'  * center: A vector of values for centering each column in X
#'  * scale: A vector of values for scaling each column in X 
#' @export
#' 
#' @examples
#' \dontrun{
#' pen <- get_data(path = "inst/extdata/penncath_lite")
#' }
#' 
#' @details
#' The .rds object should have an 'std_X' element - this is what will be used as the design matrix for analysis. This design matrix should *not* include an intercept column (this will be added later in `plmm_fit`()).
#' 
#' 
get_data <- function(path, row_id = NULL, fbm, trace = TRUE){
  
  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods 
  
  if(file.exists(bk)){
    obj <- bigsnpr::snp_attach(rdsfile = rds)
  } else {
    obj <- readRDS(rds)
  }
  
  # return data in a tractable format 
  if (missing(fbm)) {
    if (utils::object.size(obj$std_X) > 1e8) {
      warning("\nDue to the large size of X (>100 Mb), X has been returned as a file-backed matrix (FBM).
              \nTo turn this message off, explicitly specify fbm=TRUE or fbm=FALSE).")
      fbm <- TRUE
    } else {
      # if it fits, it ships 
      fbm <- FALSE
    }
  }
  
  if(!fbm){
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
    
    cat("\nReminder: the X that is returned here is column-standardized.
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
    #     \n At this time, plmm() cannot analyze design matrix X in this FBM format. Allowing such an option will
    #     \n require writing the 'meat and potatoes' of plmm() in C++, which is a work in progress. For now, 
    #     \n functions from package bigsnpr may be used for analyzing FBM data.")
    
    # std_X <- obj$std_X
    
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