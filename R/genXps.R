
#' Generate an X matrix with population structure
#'
#' This function allows you to simulate structured genetic data (SNP).
#' @param n Number of observations/samples to simulate.
#' @param nJ Number of observations in each subpopulation. The length of nJ corresponds to the number of subpopulations.
#' @param p Number of SNPs to simulate.
#' @param structureX Type of structure to simulate.
#' @param Fst The desired final FST of the admixed individuals. Ranges from 0 to 1. A high Fst implies greater differentiation among populations. Defaults to 0.1 if structureX = 1d_linear and 0.2 if structureX = indep_subpops. Otherwise defaults to NULL.
#' @param inbr Indicates whether the desired inbreeding is homogeneous or heterogeneous. Defaults to heterogeneous.
#' @param standardizeX Should the generated X matrix be standardized? Defaults to TRUE.
#' @param plot Should a plot of the kinship matrix be generated? Defaults to FALSE.
#' @param structureX_other If \code{structureX == "other"}, an matrix or SnpMatrix object with subjects in rows and SNPs in columns to be used to generate pseudophenotypes must be supplied here.
#' @param sampleCols A logical flag for whether the columns of the resultant X matrix should be scrambled. This may be desirable if the causal SNPs should change from one simulation to the next. Defaults to TRUE.
#' @export
#' 
#' @examples 
#' sim_X <- genXps(n = 200, nJ = rep(100, 2), p = 1000, structureX = "indep_subpops", inbr = "heterogeneous")

genXps <- function(n, nJ, p,
                   structureX = c("admixture", "indep_subpops", "1d_linear", "1d_circular", "independent", "other"),
                   Fst = NULL,
                   inbr = c("homogeneous", "heterogeneous"),
                   standardizeX = TRUE, plot = FALSE, structureX_other = NULL, sampleCols = TRUE){
  structureX <- match.arg(structureX)
  if (structureX == "other" & is.null(structureX_other)) stop("A matrix or SnpMatrix object must be supplied to the argument `structureX_other` if structureX == `other`")

  if (structureX == "other"){
    if (class(structureX_other)[1] == 'SnpMatrix'){
      X <- methods::as(structureX_other, "numeric")
    } else {
      X <- as.matrix(structureX_other)
    }
  } else if (structureX == "admixture"){
    dat <- utils::read.delim("https://s3.amazonaws.com/pbreheny-data-sets/admixture.txt")
    XX <- as.matrix(dat[,-1])
    polymorphic <- apply(XX, 2, stats::sd) != 0
    X <- XX[, polymorphic]
  } else if (structureX == "independent"){
    pi <- 0.2
    X <- matrix(stats::rbinom(n*p, 2, pi), nrow = n, ncol = p)
  } else if (structureX == "indep_subpops"){
    subpops <- rep(c(paste0("S", 1:length(nJ))), times = nJ)
    admix_proportions <- bnpsd::admix_prop_indep_subpops(subpops)
    if (is.null(Fst)){
      Fst <- 0.2
    }
    if (inbr == "heterogeneous"){
      inbr_subpops <- 1:length(nJ)
    } else {
      inbr_subpops <- rep(1, length(nJ))
    }
    inbr_subpops <- inbr_subpops / popkin::fst(inbr_subpops) * Fst
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    weights <- popkin::weights_subpops(subpops) # function from `popkin` package
    X <- t(bnpsd::draw_all_admix(admix_proportions = admix_proportions,
                          inbr_subpops = inbr_subpops,
                          m_loci = p)$X)
    ### plot things
    if (plot == TRUE) {
      popkin::plot_popkin(coancestry, mar = c(1, 1, 1, 3))
    }

  } else if (structureX == "1d_linear"){
    subpops <- rep(c(paste0("S", 1:length(nJ))), times = nJ)
    if (inbr == "heterogeneous"){
      inbr_subpops <- 1:length(nJ)
    } else {
      inbr_subpops <- rep(1, length(nJ))
    }
    if (is.null(Fst)){
      Fst <- 0.1
    }
    # making this smaller makes blocks more distinct
    # bigger means more blurred together and indistinguishable
    # 0.35 - 0.5 seems reasonable
    if (length(nJ) == 4){
      bias_coeff <- 0.5
    } else if (length(nJ) == 10){
      bias_coeff <- 0.5
    } else if (length(nJ) == 50){
      bias_coeff <- 0.1
    } else {
      bias_coeff <- 0.5
    }
    tmp <- bnpsd::admix_prop_1d_linear(n_ind = n,
                                k_subpops = length(nJ),
                                bias_coeff = bias_coeff,
                                coanc_subpops = inbr_subpops,
                                fst = Fst)
    admix_proportions <- tmp$admix_proportions
    inbr_subpops <- tmp$coanc_subpops
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    X <- t(bnpsd::draw_all_admix(admix_proportions = admix_proportions,
                          inbr_subpops = inbr_subpops,
                          m_loci = p)$X)
    ### plot things
    if (plot == TRUE) {
      popkin::plot_popkin(coancestry, mar = c(1, 1, 1, 3))
    }
  } else if (structureX == "1d_circular"){
    subpops <- rep(c(paste0("S", 1:length(nJ))), times = nJ)
    if (inbr == "heterogeneous"){
      inbr_subpops <- 1:length(nJ)
    } else {
      inbr_subpops <- rep(1, length(nJ))
    }
    if (is.null(Fst)){
      Fst <- 0.1
    }
    # making this smaller makes blocks more distinct
    # bigger means more blurred together and indistinguishable
    # 0.35 - 0.5 seems reasonable
    if (length(nJ) == 4){
      bias_coeff <- 0.5
    } else if (length(nJ) == 10){
      bias_coeff <- 0.5
    } else if (length(nJ) == 50){
      bias_coeff <- 0.1
    } else {
      bias_coeff <- 0.5
    }

    tmp <- bnpsd::admix_prop_1d_circular(n_ind = n,
                                       k_subpops = length(nJ),
                                       bias_coeff = bias_coeff,
                                       coanc_subpops = inbr_subpops,
                                       fst = Fst)
    admix_proportions <- tmp$admix_proportions
    inbr_subpops <- tmp$coanc_subpops
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    X <- t(bnpsd::draw_all_admix(admix_proportions = admix_proportions,
                                 inbr_subpops = inbr_subpops,
                                 m_loci = p)$X)
    ### plot things
    if (plot == TRUE) {
      popkin::plot_popkin(coancestry, mar = c(1, 1, 1, 3))
    }
  }

  # randomize column order so causal SNPs will change
  if (sampleCols) X <- X[, sample(1:ncol(X), ncol(X), replace = FALSE), drop = FALSE]
  # make numeric, not integer
  X <- apply(X, 2, as.numeric)
  # standardize
  if (standardizeX == TRUE) X <- ncvreg::std(X)
  return(X)
}
