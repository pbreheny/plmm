#' Preprocess PLINK files
#'
#' This function allows you to preprocess PLINK bed/bim/fam files for use with \code{penalizedLMM} functions. Unreliable SNPs are removed and missing values are imptued using either \code{snpStats}, or if not tagged, the HWE mean value.
#' @param prefix Character argument that is the prefix of your bed/bim/fam files.
#' @param dataDir Directory where plink files (and .sexcheck files are located if \code{sexcheck = TRUE}) are located
#' @param sexcheck Logical flag for whether or not PLINK sexcheck files should be incorporated. Defaults to FALSE. If TRUE, sexcheck files must be of the form "prefix.sexcheck"
#' @param na.strings For \code{snpStats}. Strings in .bam and .fam files to be recoded as NA. Defaults to "-9"
#' @param impute Logical flag for whether imputation should be performed. Defaults to TRUE since plmm cannot handle missing values.
#' @return A three element list object:
#' * `genotypes` The filtered and imputed genotypes in a snpMatrix object with subjects in rows and SNPs in columns.
#' * `map` A matrix of SNP data.
#' * `fam` A matrix of subject data.
#' @export
#' 
#' @examples 
#' test <- process_plink(prefix = "cad", dataDir = plink_example(path="cad.fam", parent=T))
#' 

process_plink <- function(prefix, dataDir, sexcheck = FALSE, na.strings = "-9", impute = TRUE){

  cat("\nPreprocessing", prefix, "data:\n")

  obj <- snpStats::read.plink(file.path(dataDir, prefix), na.strings = na.strings)

  # Only keep polygenic SNPs and those on chr1:22
  genotypes <- obj$genotypes[, !(snpStats::col.summary(obj$genotypes)$P.AA %in% c(NA, 1)) &
                               !(snpStats::col.summary(obj$genotypes)$P.AB %in% c(NA, 1)) &
                               !(snpStats::col.summary(obj$genotypes)$P.BB %in% c(NA, 1)) &
                               obj$map$chromosome %in% 1:22]

  cat("\nRemoving ", dim(obj$geno)[2] - dim(genotypes)[2], "SNPs that are monomorphic or outside of chromosomes 1-22.\n")

  # if sexcheck = TRUE, remove subjects with sex discrepancies
  if (sexcheck == TRUE) {
    sexdat <- data.table::fread(file.path(dataDir, paste0(prefix, ".sexcheck")))
    iid_sex_discrep <- as.character(sexdat[sexdat$STATUS == "PROBLEM", "IID"])
    if (length(iid_sex_discrep) > 0){
      genotypes <- genotypes[-match(iid_sex_discrep, rownames(genotypes)), ]
    }
    cat("\nRemoving", ifelse(length(iid_sex_discrep) > 0, length(iid_sex_discrep), 0), "subjects with sex discrepancies.\n")

  }

  if (impute) {

    # Now how many SNPs have missing data?
    missing <- table(snpStats::col.summary(genotypes)$Call.rate == 1)
    cat("\nThere are", ifelse(length(missing) == 2, missing[1], 0), "SNPs with missing data that we will attempt to impute.\n")

    # Impute missing values
    rules <- snpStats::snp.imputation(genotypes, minA=0)
    out <- snpStats::impute.snps(rules, genotypes, as.numeric=FALSE)

    # How many SNPs have missing data after imputation?
    missing <- table(snpStats::col.summary(out)$Call.rate == 1)
    cat("\nThere are", ifelse(length(missing) == 2, missing[1], 0), "SNPs with missing data after imputation.\n")

    # How many SNPs cannot be imputed and still have >= 50% missing values?
    missing <- table(snpStats::col.summary(out)$Call.rate <= 0.5)
    cat("\nOf these SNPs,", ifelse(length(missing) == 2, missing[2], 0), "have call rates <= 50% and will be removed.\n")

    # Throw out SNPs that have >= 50% missingness, even after imputation
    out2 <- out[, snpStats::col.summary(out)$Call.rate > 0.5]

    # How many SNPs will be imputed with mean values?
    missing <- table(snpStats::col.summary(out2)$Call.rate == 1)
    cat("\nThere are", ifelse(length(missing) == 2, missing[1], 0), "SNPs that still have missing data which will be imputed with the HWE expected value.\n")

    # keep snp order for later
    order_snps <- colnames(out2)

    # which SNPs have missingness
    to_impute <- which(snpStats::col.summary(out2)$Call.rate < 1)
    miss <- out2[, to_impute]

    imputed_mean <- apply(methods::as(miss, "numeric"), 2, function(s){
      these <- which(is.na(s))
      s[these] <- mean(s, na.rm = TRUE)
      s <- snpStats::mean2g(s)
      return(s)
    })

    out3 <- out2
    out3@.Data[, to_impute] <- imputed_mean
    genotypes <- out3
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

  # keep corresponding map data
  map <- obj$map[colnames(genotypes),]

  # keep corresponding fam data
  fam <- obj$fam[rownames(genotypes),]
  fam$prefix <- prefix

  # done...
  cat("\nPreprocessing", prefix, "data DONE!\n",
      "\nSubjects:", nrow(genotypes),
      "\nSNPs:", ncol(genotypes),
      "\nMissing values:", missing, "\n")

  return(list(genotypes = genotypes,
              map = map,
              fam = fam))
}

