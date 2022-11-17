#' Preprocess PLINK files
#'
#' This function allows you to preprocess PLINK bed/bim/fam files for use with \code{penalizedLMM} functions. Unreliable SNPs are removed and missing values are imptued using either \code{snpStats}, or if not tagged, the HWE mean value.
#' @param prefix Character argument that is the prefix of your bed/bim/fam files.
#' @param dataDir Directory where plink files (and .sexcheck files are located if \code{sexcheck = TRUE}) are located
#' @param sexcheck Logical flag for whether or not PLINK sexcheck files should be incorporated. Defaults to FALSE. If TRUE, sexcheck files must be of the form "prefix.sexcheck"
#' @param na.strings For \code{snpStats}. Strings in .bam and .fam files to be recoded as NA. Defaults to "-9"
#' @param impute Logical flag for whether imputation should be performed. Defaults to TRUE since plmm cannot handle missing values.
#' @param quiet Logical flag: should progress messages be printed? Defaults to FALSE. 
#' @return A three element list object:
#' * `genotypes` The filtered and imputed genotypes in a snpMatrix object with subjects in rows and SNPs in columns. If coerce=T, this is a numeric matrix. 
#' * `map` A matrix of SNP data.
#' * `fam` A matrix of subject data.
#' @export
#' 
#' @examples 
#' test <- process_plink(prefix = "cad_lite", dataDir = plink_example(path="cad_lite.fam", parent=T))
#' \dontrun{
#' test2 <- process_plink(prefix = "cad", dataDir = plink_example(path="cad.fam", parent=T))
#' }
#' 
#' 

process_plink <- function(prefix, dataDir, sexcheck = FALSE, na.strings = "-9", impute = TRUE, quiet=FALSE){

  if(!quiet){
    cat("\nPreprocessing", prefix, "data:\n")
  }

  obj <- snpStats::read.plink(file.path(dataDir, prefix), na.strings = na.strings)
  
  # save snpStats::colsummary() result - this is a data frame with numeric columns
  obj_col_summary <- snpStats::col.summary(obj$genotypes)
  
  # coerce types for 'genotypes'
  genotypes <- as(obj$genotypes, 'numeric') # coerce SnpMatrix into a numeric matrix 
  
  
  # only consider SNPs on chromosomes 1-22
  map <- subset(x = obj$map, subset = chromosome %in% 1:22)
  genotypes <- genotypes[,rownames(map)]
  obj_col_summary <- obj_col_summary[colnames(genotypes),]
  
  # only keep polygenic SNPs 
  obj_col_summary <- subset(x = obj_col_summary,
                            # step 1: remove SNPs that are missing a genotype
                            subset = !(is.na(P.AA) | is.na(P.AB) | is.na(P.BB)))
  
  obj_col_summary <- subset(x = obj_col_summary,
                            # step 2: remove SNPs that have no variation (i.e. only one genotype is present)
                            subset = !(P.AA == 1 | P.AB == 1 | P.BB == 1))
  
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

