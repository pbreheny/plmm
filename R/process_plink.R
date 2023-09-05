#' Preprocess PLINK files using the `bigsnpr` package
#' 
#' @param data_dir The path to the bed/bim/fam data files 
#' @param prefix The prefix (as a character string) of the bed/fam data files 
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param impute_method If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are: 
#'  * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details. 
#'  * random: Imputes sampling according to allele frequencies.
#'  * mean0: Imputes the rounded mean.
#'  * mean2: Imputes the mean rounded to 2 decimal places.
#'  * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set. 
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param ... Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' process_plink(data_dir = "../temp_files",
#'  prefix = "penncath_lite",
#'   impute_method = "xgboost",
#'    quiet = F)
#' }
process_plink <- function(data_dir,
                          prefix,
                          impute = TRUE,
                          impute_method = c('mode'),
                          quiet = FALSE,
                          gz = FALSE,
                          outfile,
                          ...){
  
  # start log ------------------------------------------
  if(missing(outfile)){
    outfile = "process_plink.log"
    } else {
      outfile = paste0(outfile, ".log")
    }
  log_con <- file(outfile)
  cat("### Processing PLINK files for PLMM ###", file = log_con)
  cat("\nLogging to ", outfile, file = outfile, append = TRUE)
  cat("\nPreprocessing", prefix, "data:", file = outfile, append = TRUE)
  
  if(!quiet){
    cat("\nLogging to", outfile)
    cat("\nPreprocessing", prefix, "data:")
  }
  
  # read in PLINK files --------------------------------
  path <- paste0(data_dir, "/", prefix, ".rds")
  
  
  # Create the RDS file first 
  cat("\nCreating ", prefix, ".rds\n", file = outfile, append = TRUE)
  if(!quiet){
    cat("\nCreating ", prefix, ".rds\n")
    
    # check for compressed files 
    if (gz){
      cat("\nUnzipping .gz files - this could take a second", file = outfile, append = TRUE)
      if (!quiet){cat("\nUnzipping .gz files - this could take a second")}
      system(paste0("gunzip -k ", file.path(data_dir, paste0(prefix, "*"))))
    }
    
    bigsnpr::snp_readBed(bedfile = paste0(data_dir, "/", prefix, ".bed"))
    obj <- bigsnpr::snp_attach(path)
  }
  
  # chromosome check ---------------------------------
  # only consider SNPs on chromosomes 1-22
  chr_range <- range(obj$map$chromosome)
  if(chr_range[1] < 1 | chr_range[2] > 22){
    cat("PLMM only analyzes autosomes -- removing chromosomes outside 1-22")
    cat("PLMM only analyzes autosomes -- removing chromosomes outside 1-22",
        file = outfile, append = TRUE)
    
    original_dim <- dim(obj$genotypes)[2]
    chr_filtered <- bigsnpr::snp_subset(obj,
                                        ind.col = obj$map$chromosome %in% 1:22)
    obj <- bigsnpr::snp_attach(chr_filtered)
    new_dim <- dim(obj$genotypes)[2]
    
    cat("\nRemoved ", original_dim - new_dim, "SNPs that are outside of chromosomes 1-22.",
        file = outfile, append = TRUE)
    if(!quiet){
      cat("\nRemoved ", original_dim - new_dim, "SNPs that are outside of chromosomes 1-22.")
      
    }
  }
  
  # TODO: figure out how to add a 'sexcheck' with bigsnpr functions
  # e.g., if sexcheck = TRUE, remove subjects with sex discrepancies
   
  
  chr <- obj$map$chromosome
  X   <- obj$genotypes
  pos <- obj$map$physical.pos
  
  
  # save these counts (like 'col_summary' obj from snpStats package)
  counts <- bigstatsr::big_counts(X) # NB this is a matrix 
  
  # identify monomorphic SNPs --------------------------------
  constants_idx <- apply(X = counts[1:3,],
                         MARGIN = 2,
                         # see which ~called~ features have all same value
                         FUN = function(c){sum(c == sum(c)) > 0})
  
  cat("\nThere are ", sum(constants_idx), " constant features in the data",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are ", sum(constants_idx), " constant features in the data")
  }
  
  # notify about missing values ----------------------------
  na_idx <- counts[4,] > 0
  prop_na <- counts[4,]/nrow(X)
  
  cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values",
      file = outfile, append = TRUE)
  cat("\nOf these, ", sum(prop_na > 0.5),
      " are missing in at least 50% of the samples",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values")
    cat("\nOf these, ", sum(prop_na > 0.5), " are missing in at least 50% of the samples")
  }
  
  # imputation ------------------------------------------
  if(!quiet & impute){
    # catch for misspellings
    if(!(impute_method %in% c('mode', 'random', 'mean0', 'mean2', 'xgboost'))){
      stop("\nImpute method is misspecified or misspelled. Please use one of the 
           \n5 options listed in the documentation.")
    }
    cat("\nImputing the missing values using ", impute_method, " method\n")
  }
  
  if(impute){
    cat("\nImputing the missing values using ", impute_method, " method",
        file = outfile, append = TRUE)
    
    if(impute_method %in% c('mode', 'random', 'mean0', 'mean2')){ 
       # NB: this will overwrite obj$genotypes
    obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                   ncores = bigstatsr::nb_cores(),
                                                   method = impute_method) # dots can pass other args
      } else if (impute_method == "xgboost"){

      imp <- bigsnpr::snp_fastImpute(Gna = X,
                                     ncores = bigstatsr::nb_cores(),
                                     infos.chr = chr,
                                     seed = as.numeric(Sys.Date()),
                                     ...) # dots can pass other args
      
      cat("\n ***************** NOTE ********************************
          \n August 2023: With the xgboost imputation method, there have been some issues (particularly
          \n on Mac OS) with warnings that appear saying 'NA or NaN values in the 
          \n resulting correlation matrix.' However, we (penalizedLMM authors) have
          \n not seen missing values appear in the results -- the imputed data
          \n does not show any NA or NaN values, and models fit on these data run without issue. 
          \n We are actively investigating this warning message, and will
          \n make a note in a future release. If using xgboost, proceed with 
          \n caution and file an issue if you notice any problems downstream.
          \n ********************************************************")

      # save imputed values (NB: will overwrite obj$genotypes)
      obj$genotypes$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }
    # now, save the new object -- this will have imputed values and constants_idx
    obj$constants_idx <- constants_idx
    obj <- bigsnpr::snp_save(obj)
    
    cat("\nDone with imputation. File formatting in progress.",
        file = outfile, append = TRUE)
    
  }
  
  if(!quiet & impute){cat("\nDone with imputation. Processed files now saved as .rds object.")}
  close(log_con)
}



