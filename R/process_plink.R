#' Preprocess PLINK files using the `bigsnpr` package
#' 
#' @param data_dir The path to the bed/bim/fam data files 
#' @param prefix The prefix (as a character string) of the bed/fam data files 
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param impute_method If 'impute' = TRUE, this argument will specify the kind of imputation desired. This is passed directly to `bigsnpr::snp_fastImputeSimple()`. Defaults to 'mode'. 
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' #' Note: a more nuanced 'impute' argument is still under construction. Only "simple" method is available at this time, but we hope to add "xgboost" as an option in the future.
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' lite <- process_plink(data_dir = plink_example(path = "penncath_lite.bed.gz", parent = T), prefix = "penncath_lite", gz = TRUE)
#' mid <- process_plink(data_dir = "/Users/tabithapeter/Desktop/penalizedLMM/data-raw", prefix = "penncath_mid", gz = F, method = "mean2")
#' }
process_plink <- function(data_dir,
                          prefix,
                          impute = TRUE,
                          impute_method = 'mode',
                          quiet = FALSE,
                          gz = FALSE,
                          outfile){
  
  # start log 
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
  
  # read in PLINK files 
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
  
  # identify monomorphic SNPs 
  constants_idx <- apply(X = counts[1:3,],
                         MARGIN = 2,
                         # see which ~called~ features have all same value
                         FUN = function(c){sum(c == sum(c)) > 0})
  
  cat("\nThere are ", sum(constants_idx), " constant features in the data",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are ", sum(constants_idx), " constant features in the data")
  }
  
  # notify about missing values
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
  
  if(!quiet & impute){
    cat("\nImputing the missing values using ", impute_method, " method\n")
  }
  
  if(impute){
    cat("\nImputing the missing values using ", impute_method, " method",
        file = outfile, append = TRUE)
    
    # NB: this will overwrite obj$genotypes
    obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                   ncores = bigstatsr::nb_cores(),
                                                   method = impute_method) # dots can pass other args
    
    
    # TODO: come back here and try to get the 'xgboost' method to work
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
    
    # now, save the new object -- this will have imputed values and constants_idx
    obj$constants_idx <- constants_idx
    obj <- bigsnpr::snp_save(obj)
    
    cat("\nDone with imputation. File formatting in progress.",
        file = outfile, append = TRUE)
    
  }
  
  if(!quiet & impute){cat("\nDone with imputation. Processed files now saved as .rds object.")}
  close(log_con)
}



#' This function allows you to preprocess PLINK bed/bim/fam files using `snpStats`, an older package compared to `bigsnpr`. Unreliable SNPs are removed and missing values are imptued using either \code{snpStats}, or if not tagged, the HWE mean value.
#' @param prefix Character argument that is the prefix of your bed/bim/fam files.
#' @param dataDir Directory where plink files (and .sexcheck files are located if \code{sexcheck = TRUE}) are located
#' @param gz A logical: are the PLINK files g-zipped (i.e. are the file endings '.gz')?. Defaults to FALSE. NOTE: if TRUE, process_plink_snpStats will unzip your zipped files 
#' @param sexcheck Logical flag for whether or not PLINK sexcheck files should be incorporated. Defaults to FALSE. If TRUE, sexcheck files must be of the form "prefix.sexcheck" in the same directory as the bed/fam files
#' @param na.strings For \code{snpStats}. Strings in .bim and .fam files to be recoded as NA. Defaults to "-9"
#' @param impute Logical flag for whether imputation should be performed. Defaults to TRUE since plmm cannot handle missing values.
#' @param quiet Logical flag: should progress messages be printed? Defaults to FALSE. 
#' @return A three element list object:
#' * `genotypes` The filtered and imputed genotypes in a snpMatrix object with subjects in rows and SNPs in columns. This is a numeric matrix. 
#' * `map` A matrix of SNP data.
#' * `fam` A matrix of subject data.
#' 
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' # the output of calls to 'plink_example' will vary according to the users' directory structure
#' test <- process_plink_snpStats(prefix = "cad_lite", dataDir = plink_example(path="cad_lite.fam.gz", parent=T), gz = TRUE)
#' test2 <- process_plink_snpStats(prefix = "cad_mid", dataDir = plink_example(path="cad_mid.fam.gz", parent=T))
#' }
#' 
#' 

process_plink_snpStats <- function(prefix,
                          dataDir,
                          gz = FALSE,
                          sexcheck = FALSE,
                          na.strings = "-9",
                          impute = TRUE,
                          quiet=FALSE){
  
  if(!quiet){
    cat("\nPreprocessing", prefix, "data:\n")
  }
  
  # check for compressed files 
  if (gz){
    if (!quiet){cat("Unzipping .gz files - this could take a second")}
    system(paste0("gunzip -k ", file.path(dataDir, paste0(prefix, "*"))))
  }
  
  obj <- snpStats::read.plink(file.path(dataDir, prefix), na.strings = na.strings)
  
  # save snpStats::colsummary() result - this is a data frame with numeric columns
  obj_col_summary <- snpStats::col.summary(obj$genotypes)
  
  # coerce types for 'genotypes'
  genotypes <- as(obj$genotypes, 'numeric') # coerce SnpMatrix into a numeric matrix 
  
  
  # only consider SNPs on chromosomes 1-22
  map <- obj$map[obj$map$chromosome %in% 1:22,]
  genotypes <- genotypes[,rownames(map)]
  obj_col_summary <- obj_col_summary[colnames(genotypes),]
  
  # only keep polygenic SNPs 
  
  # step 1: remove SNPs that are missing a genotype
  missing_genotypes <- c(is.na(obj_col_summary$P.AA) |
                           is.na(obj_col_summary$P.AB) |
                           is.na(obj_col_summary$P.BB))
  obj_col_summary <- obj_col_summary[!missing_genotypes,]
  # step 2: remove SNPs that have no variation (i.e. only one genotype is present)
  no_var <- c(obj_col_summary$P.AA == 1 | 
                obj_col_summary$P.AB == 1 | 
                obj_col_summary$P.BB == 1)
  obj_col_summary <- obj_col_summary[!no_var,]
  
  # subset the genotypes data to remove the monogenic SNPs
  genotypes <- genotypes[, rownames(obj_col_summary)]
  
  if(!quiet){
    cat("\nRemoved ", dim(obj$genotypes)[2] - dim(genotypes)[2], "SNPs that are monomorphic or outside of chromosomes 1-22.\n")
    
  }
  
  # if sexcheck = TRUE, remove subjects with sex discrepancies
  if (sexcheck == TRUE) {
    sexdat <- data.table::fread(file.path(dataDir, paste0(prefix, ".sexcheck")))
    iid_sex_discrep <- as.character(sexdat[sexdat$STATUS == "PROBLEM", "IID"])
    if (length(iid_sex_discrep) > 0){
      genotypes <- genotypes[-match(iid_sex_discrep, rownames(genotypes)), ]
    }
    if(!quiet){
      cat("\nRemoving", ifelse(length(iid_sex_discrep) > 0, length(iid_sex_discrep), 0), "subjects with sex discrepancies.\n")
    }
    
  }
  
  if (impute) {
    
    # Now how many SNPs have missing data?
    missing <- table(obj_col_summary$Call.rate == 1)
    if(!quiet){
      cat("\nThere are", ifelse(length(missing) == 2, missing[1], 0), "SNPs with missing data that we will attempt to impute.\n")
      
    }
    
    # Impute missing values
    rules <- snpStats::snp.imputation(obj$genotypes[,colnames(genotypes)], minA=0)
    Imp <- snpStats::impute.snps(rules = rules, snps = obj$genotypes[,colnames(genotypes)])
    to_fill_in <- is.na(genotypes) # these are the indices marking which values to 'fill in' with the imputation 
    genotypes[to_fill_in] <- Imp[to_fill_in] # this step fills in the missing values with imputed values
    
    # How many SNPs have missing data after imputation?
    still_na <- apply(genotypes, MARGIN = 2, FUN = function(x){sum(is.na(x))})
    if(!quiet){
      cat("\nThere are", sum(still_na > 0), "SNPs with missing data after imputation.\n")
      
    }
    
    # define call rates
    call_rates <- obj_col_summary$Call.rate
    
    # how many SNPs cannot be imputed and still have >= 50% missing values?
    above_half_na <- sum(call_rates < 0.5)
    if(!quiet){
      cat("\nOf these SNPs,", above_half_na, "have call rates <= 50% and will be removed.\n")
      
    }
    
    # Throw out SNPs that have >= 50% missingness, even after imputation
    genotypes_keep <- genotypes[, which(call_rates > 0.5)]
    
    # How many SNPs will be imputed with mean values?
    imp_hwe <- apply(genotypes_keep, MARGIN = 2, FUN = function(x){sum(is.na(x))})
    if(!quiet){
      cat("\nThere are", sum(imp_hwe > 0), "SNPs that still have missing data which will be imputed with the HWE expected value.\n")
      
    }
    
    # keep snp order for later
    order_snps <- colnames(genotypes_keep)
    
    # which SNPs have missingness
    call_rates_keep <- apply(genotypes_keep, MARGIN = 2, FUN = function(x){sum(!is.na(x))/length(x)})
    to_impute <- which(call_rates_keep < 1)
    Imp_hwe <- genotypes_keep[, to_impute]
    
    imputed_mean <- apply(Imp_hwe, 2, function(s){
      na_vals <- which(is.na(s))
      s[na_vals] <- mean(s, na.rm = TRUE)
      
      return(s)
    })
    
    genotypes_none_missing <- genotypes_keep
    genotypes_none_missing[, to_impute] <- imputed_mean
    genotypes <- genotypes_none_missing
  }
  
  # make sure there are no missing values remaining / count missing values
  missing <- 0
  chunks <- ceiling(nrow(genotypes) / 100)
  start <- 1
  for (i in 1:chunks){
    stop <- min(i*100, nrow(genotypes))
    missing <- missing + sum(is.na(genotypes[start:stop]))
    start <- stop + 1
  }
  
  # triple check ordering: 
  if(all.equal(order_snps, colnames(genotypes))){cat("Imputation complete; alignment check passed")} else {
    warning("Alignment issue has occured - SNPs in supplied data and imputed data do not match.")
  }
  
  # keep corresponding map data
  map <- obj$map[colnames(genotypes),]
  
  # keep corresponding fam data
  fam <- obj$fam[rownames(genotypes),]
  fam$prefix <- prefix
  
  
  # done...
  if(!quiet){
    cat("\nPreprocessing", prefix, "data DONE!\n",
        "\nSubjects:", nrow(genotypes),
        "\nSNPs:", ncol(genotypes),
        "\nMissing values:", missing, "\n")
  }
  
  
  return(list(genotypes = genotypes,
              map = map,
              fam = fam))
}



