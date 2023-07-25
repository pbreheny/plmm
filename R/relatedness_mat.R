#' Calculate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix. If a numeric matrix is supplied, the RRM (Hayes, 2009) is used
#' and is computed XX'/p, where X is standardized. 
#' @param X Either a numeric matrix of genotypes (subjects in rows, SNPs in columns), or a character prefix for PLINK bed/bim/fam files.
#' @param ... Optional arguments in the case that a character prefix is specified.
#' @importFrom SNPRelate snpgdsBED2GDS snpgdsGRM snpgdsOpen snpgdsClose
#' @export
#' 
#' @examples 
#' RRM <- relatedness_mat(X = admix$X)

relatedness_mat <- function(X, ...){
  UseMethod("relatedness_mat")
}


#' @export
relatedness_mat.character <- function(X, ...){
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
relatedness_mat.matrix <- function(X, ...){
  # NB: X is standardized as part of the RRM calculation
  rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
  return(rrm)
}
