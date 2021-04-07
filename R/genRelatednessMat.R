#' Generate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix. If a numeric matrix is supplied, the RRM (Hayes, 2009) is used
#' and is computed XX'/p. If a character argument which describes the location and prefix of PLINK bed/bim/bam files is supplied,
#' a GRM is computed using the GCTA method of SNPrelate, unless another method is specified as an optional argument.
#' @param X Either a numeric matrix of genotypes (subjects in rows, SNPs in columns), or a character prefix for PLINK bed/bim/fam files.
#' @param ... Optional arguments in the case that a character prefix is specified.
#' @importFrom SNPRelate snpgdsBED2GDS snpgdsGRM snpgdsOpen snpgdsClose
#' @export
genRelatednessMat <- function(X, ...){
  UseMethod("genRelatednessMat")
}


#' @export
genRelatednessMat.character <- function(X, ...){
  args <- list(...)
  gds.fn <- lapply(c(bed='bed', bim='bim', fam='fam', gds='gds'), function(y) paste0(X, '.', y))
  snpgdsBED2GDS(gds.fn$bed, gds.fn$fam, gds.fn$bim, gds.fn$gds)
  genofile <- snpgdsOpen(gds.fn$gds)
  if (!('method' %in% names(args))) args$method <- 'GCTA'
  args$gdsobj <- genofile
  grm <- do.call('snpgdsGRM', args)
  snpgdsClose(genofile)
  unlink(paste0(X, '.gds'), force=TRUE)
  return(grm)
}


#' @export
genRelatednessMat.matrix <- function(X, ...){
  rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
  return(rrm)
}
