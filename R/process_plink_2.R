#' An updated function for processing PLINK files
#' 
#' @param data_dir The path to the bed/bim/fam data files 
#' @param prefix The prefix (as a character string) of the bed/fam data files 
#' @param rds Logical: does an rds file for these data files already exist in \code{data_dir}? Defaults to FALSE
#' @param impute Character string indicating type of imputation to use. "simple" will use \code{bigsnpr::snp_fastImputeSimple}, and "xgboost" will use \code{bigsnpr::snp_fastImpute}. Defaults to "simple."
#' Note: 'impute' argument is still under construction. Only "simple" method is available at this time.
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param ... Other arguments to bigsnpr::snp_fastImpute or bigsnpr::snp_fastImputeSimple (depending on choice of \code{impute})
#' 
#' @examples 
#' 
#' \dontrun{
#' cad_mid <- process_plink_2(data_dir = plink_example(path = "cad_mid.bed", parent = T), prefix = "cad_mid", rds = T, impute = "simple", method = "mean0")
#' }
process_plink_2 <- function(data_dir, prefix, rds = FALSE, impute = "simple", quiet = FALSE, ...){
  
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
      cat("\nRemoved ", old_dim - new_dim, "SNPs that are outside of chromosomes 1-22.\n")
      
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
    cat("\nImputing the missing values using ", impute, "method\n")
  }
  
  # impute missing values
  if(impute == "simple"){
    obj$imputed_X <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                   ncores = bigstatsr::nb_cores(),
                                                   ...) # dots can pass method
    
  } else if (impute == "xgboost"){
    # FIXME: this is still not working; when I kill it, I get this warning: 
    # "NA or NaN values in the resulting correlation matrix."
    obj$imputed_X <- bigsnpr::snp_fastImpute(Gna = X,
                                             ncores = bigstatsr::nb_cores(),
                                             infos.chr = chr,
                                             seed = as.numeric(Sys.Date()),
                                             ...) # dots can pass method
    
    
  } else stop("Argument impute must be either simple or xgboost.")
  
  # now, to make the save the imputed values
  obj$imputed_X$code256 <- bigsnpr::CODE_IMPUTE_PRED
  obj <- bigsnpr::snp_save(obj)
  
  # return data in a tractable format 
  # NB: this assumes that X will fit into memory!! 
  # TODO: add nuance here 
  imputed_X <- obj$imputed_X[,]
  
  if(!quiet){cat("\nDone!\n")}
  
  return(list(X = imputed_X,
              constants_idx = which(constants_idx == "TRUE")))
  
  
}

