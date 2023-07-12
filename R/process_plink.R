#' Preprocess PLINK files using the `bigsnpr` package
#' 
#' @param data_dir The path to the bed/bim/fam data files 
#' @param prefix The prefix (as a character string) of the bed/fam data files 
#' @param rds Logical: does an rds file for these data files already exist in \code{data_dir}? Defaults to FALSE
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param row_id Character string indicating which IDs to use for the rownames of the genotype matrix. Can choose "fid" or "iid", corresponding to the first or second columns in the PLINK .fam file. Defaults to NULL. 
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param ... Other arguments to bigsnpr::snp_fastImpute or bigsnpr::snp_fastImputeSimple (depending on choice of \code{impute})
#' #' Note: a more nuanced 'impute' argument is still under construction. Only "simple" method is available at this time, but we hope to add "xgboost" as an option in the future.
#' 
#' @return A list of two components: 
#' * X: the fully-imputed design matrix, whose columns are the features and whose rows are the observations. 
#' * constants_idx: A numeric vector with the indices of the constant features. An index of 5 indicates that the 5th feature is constant (i.e. has zero variance)
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' lite <- process_plink(data_dir = plink_example(path = "penncath_lite.bed.gz", parent = T), prefix = "penncath_lite", rds = F, gz = TRUE)
#' mid <- process_plink(data_dir = "/Users/tabithapeter/Desktop/penalizedLMM/data-raw", prefix = "penncath_mid", rds = F, gz = F, row_id = "fid", method = "mode")
#' }
process_plink <- function(data_dir,
                            prefix,
                            rds = FALSE,
                            impute = TRUE,
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
  if(impute){
    # NB: this will overwrite obj$genotypes
    obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                   ncores = bigstatsr::nb_cores(),
                                                   ...) # dots can pass method (mean, mode, etc.)
    
    # now, save the imputed values
    obj <- bigsnpr::snp_save(obj)
  }
  
  
  
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



#' This function allows you to preprocess PLINK bed/bim/fam files using `snpStats`, an older package compared to `bigsnpr`. Unreliable SNPs are removed and missing values are imptued using either \code{snpStats}, or if not tagged, the HWE mean value.
#' @param prefix Character argument that is the prefix of your bed/bim/fam files.
#' @param dataDir Directory where plink files (and .sexcheck files are located if \code{sexcheck = TRUE}) are located
#' @param gz A logical: are the PLINK files g-zipped (i.e. are the file endings '.gz')?. Defaults to FALSE. NOTE: if TRUE, process_plink_snpStats will unzip your zipped files 
#' @param sexcheck Logical flag for whether or not PLINK sexcheck files should be incorporated. Defaults to FALSE. If TRUE, sexcheck files must be of the form "prefix.sexcheck" in the same directory as the bed/fam files
#' @param na.strings For \code{snpStats}. Strings in .bam and .fam files to be recoded as NA. Defaults to "-9"
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



