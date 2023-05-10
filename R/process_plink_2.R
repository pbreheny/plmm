#' An updated function for processing PLINK files
#' 
#' @param data_dir The path to the bed/bim/fam data files 
#' @param prefix The prefix (as a character string) of the bed/fam data files 
#' @param rds Logical: does an rds file for these data files already exist in \code{data_dir}? Defaults to FALSE
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param row_id Character string indicating which IDs to use for the rownames of the genotype matrix. Can choose "fid" or "iid", corresponding to the first or second columns in the PLINK .fam file. Defaults to NULL. 
#' @param ... Other arguments to bigsnpr::snp_fastImpute or bigsnpr::snp_fastImputeSimple (depending on choice of \code{impute})
#' #' Note: an 'impute' argument is still under construction. Only "simple" method is available at this time, but we hope to add "xgboost" as an option in the future.
#' 
#' @return A list of two components: 
#' * X: the fully-imputed design matrix, whose columns are the features and whose rows are the observations. 
#' * constants_idx: A numeric vector with the indices of the constant features. An index of 5 indicates that the 5th feature is constant (i.e. has zero variance)
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' lite <- process_plink_2(data_dir = plink_example(path = "penncath_lite.bed.gz", parent = T), prefix = "penncath_lite", rds = F, gz = TRUE)
#' mid <- process_plink_2(data_dir = "/Users/tabithapeter/Desktop/penalizedLMM/data-raw", prefix = "penncath_mid", rds = F, gz = F, row_id = "fid", method = "mode")
#' }
process_plink_2 <- function(data_dir,
                            prefix,
                            rds = FALSE,
                            impute = "simple",
                            quiet = FALSE,
                            gz = FALSE,
                            row_id = NULL,
                            ...){
  
  if(!quiet){
    cat("\nPreprocessing", prefix, "data:\n")
  }
  
  # read in PLINK files 
  path <- paste0(data_dir, "/", prefix, ".rds")
  
  if(rds){
    # if an RDS file already exists, use it 
    obj <- bigsnpr::snp_attach(path)
  } else {
    # else create the RDS file first 
    if(!quiet){
      cat("\nCreating ", prefix, ".rds\n")
    }

    # check for compressed files 
    if (gz){
      if (!quiet){cat("Unzipping .gz files - this could take a second")}
      system(paste0("gunzip -k ", file.path(data_dir, paste0(prefix, "*"))))
    }
    
    bigsnpr::snp_readBed(bedfile = paste0(data_dir, "/", prefix, ".bed"))
    obj <- bigsnpr::snp_attach(path)
  }
  
  # only consider SNPs on chromosomes 1-22
  chr_range <- range(obj$map$chromosome)
  
  if(chr_range[1] != 1 | chr_range[2] != 22){
    original_dim <- dim(obj$genotypes)[2]
    chr_filtered <- bigsnpr::snp_subset(obj, ind.col = chr %in% 1:22)
    obj <- bigsnpr::snp_attach(chr_filtered)
    new_dim <- dim(obj$genotypes)[2]
    
    if(!quiet){
      cat("\nRemoved ", original_dim - new_dim, "SNPs that are outside of chromosomes 1-22.\n")
      
    }
  }
  
  
  # if sexcheck = TRUE, remove subjects with sex discrepancies
  # TODO: figure out how to do this with bigsnpr functions 
  
  chr <- obj$map$chromosome
  X   <- obj$genotypes
  pos <- obj$map$physical.pos
  
  
  # save these counts (like 'col_summary' obj from snpStats package)
  counts <- bigstatsr::big_counts(X) # NB this is a matrix 
  
  # identify monomorphic SNPs 
  constants_idx <- apply(X = counts[1:3,],
                         MARGIN = 2,
                         # see which ~called~ features have all same value
                         FUN = function(c){sum(c == sum(c)) > 0})
  if(!quiet){
    cat("\nThere are ", sum(constants_idx), " constant features in the data\n")
  }
  
  # notify about missing values
  na_idx <- counts[4,] > 0
  prop_na <- counts[4,]/nrow(X)
  if(!quiet){
    cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values\n")
    cat("\nOf these, ", sum(prop_na > 0.5), " are missing in at least 50% of the samples\n")
  }
  
  if(!quiet){
    # TODO: add detail (mean, mode, etc.) to the output here
    # impute_type <- ifelse(missing(method), impute, method)
    cat("\nImputing the missing values using ",impute, "method\n")
  }
  
  # impute missing values
  # NB: this will overwrite obj$genotypes
  obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                 ncores = bigstatsr::nb_cores(),
                                                 ...) # dots can pass method (mean, mode, etc.)
  
  # now, save the imputed values
  obj <- bigsnpr::snp_save(obj)
  
  
  # TODO: come back here and try to get the 'xgboost' method to work
  # if(impute == "simple"){
  #   # NB: this will overwrite obj$genotypes
  #   obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
  #                                                  ncores = bigstatsr::nb_cores(),
  #                                                  ...) # dots can pass method
  #   
  #   # now, save the imputed values
  #   obj <- bigsnpr::snp_save(obj)
  # 
  #   
  # } else if (impute == "xgboost"){
  #   imp <- bigsnpr::snp_fastImpute(Gna = X,
  #                                  ncores = bigstatsr::nb_cores(),
  #                                  infos.chr = chr,
  #                                  seed = as.numeric(Sys.Date()),
  #                                  ...) # dots can pass method
  #   
  #   # save imputed values (NB: will overwrite obj$genotypes)
  #   obj$genotypes$code256 <- bigsnpr::CODE_IMPUTE_PRED
  #   obj <- bigsnpr::snp_save(obj)
  # 
  #   
  # } else stop("Argument impute must be either simple or xgboost.")
  # 
  # 
  
  
  if(!quiet){cat("\nDone!\n")}
  
  # return data in a tractable format 
  # NB: this assumes that X will fit into memory!! 
  # TODO: add nuance here 
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
              constants_idx = which(constants_idx == "TRUE")))
  
  # TODO: create a .log file? 
  
}

