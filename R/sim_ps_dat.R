
#' Simulate data with population structure
#'
#' This function allows you to simulate structured genetic data (SNP) with an unobserved environmental confounding effect.
#' @param n Number of observations/samples to simulate. Defaults to 200. 
#' @param p Number of SNPs to simulate. Defaults to 1000. 
#' @param p1 Number of SNPs that are causal. Defaults to floor(p/2).
#' @param nJ Number of observations in each subpopulation. The length of nJ corresponds to the number of subpopulations. Defaults to rep(50, 4). 
#' @param structureX Type of structure to simulate.
#' @param Fst The desired final FST of the admixed individuals. Ranges from 0 to 1. A high FST implies greater differentiation among populations. Defaults to 0.1 if structureX = 1d_linear and 0.2 if structureX = indep_subpops. Otherwise defaults to NULL.
#' @param inbr Indicates whether the desired inbreeding is homogeneous or heterogeneous. Defaults to heterogeneous.
#' @param structureGamma The desired structure of the environmental confounding effect. Defaults to 'dichotomous_discordant'
#' @param eta The desired proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR.
#' @param xi The desired proportion of the non-signal variance in the outcome that is attributable to unobserved environmental confounding effects.
#' @param standardizeX Should the generated X matrix be standardized? Defaults to TRUE.
#' @param structureX_other If \code{structureX == "other"}, an matrix or SnpMatrix object with subjects in rows and SNPs in columns to be used to generate pseudophenotypes must be supplied here.
#' @param sampleCols A logical flag for whether the columns of the resultant X matrix should be scrambled. This may be desirable if the causal SNPs should change from one simulation to the next. Defaults to TRUE.
#' @param B0 Optional. Additional intercept value.
#' 
#' @return A list with 7 elements: a matrix of SNP data (X), a single-column matrix of outcome values (y), a single column matrix 'env', a vector of coefficients (beta), a matrix of Z values allocating environmental effects among subjects (Z), a vector of numeric values representing environmental effects (gamma), and the type of structure used in the SNP data (structureX)
#' 
#' @export
#' 
#' @examples 
#' sim_dat <- sim_ps_dat(structureX = "1d_linear")
#' example_fit <- plmm(sim_dat$X, sim_dat$y, K = sim_dat$X%*%t(sim_dat$X))


sim_ps_dat <- function(n = 200, p = 1000, p1 = floor(p/2), nJ = rep(50, 4),
                      # structureX = c("admixture", "indep_subpops", "1d_linear", "1d_circular", "independent", "other"),
                      structureX = "indep_subpops",
                      Fst = NULL,
                      # inbr = c("heterogeneous", "homogeneous"),
                      inbr = "heterogeneous",
                      structureGamma = c("dichotomous_discordant"),
                      eta = 0.8,
                      xi = 0,
                      standardizeX = TRUE,
                      structureX_other = NULL,
                      sampleCols = TRUE,
                      B0){

  # structureX <- match.arg(structureX)
  if (structureX == "other" & is.null(structureX_other)) stop("A matrix or SnpMatrix object must be supplied to the argument `structureX_other` if structureX == `other`")
  J <- length(nJ)

  # dimension checks
  if (structureX == "admixture"){
    if ((n != 197 || p != 98 || !missing(nJ) & length(nJ) != 4)){
      warning("Admixture structureX is specified which will override supplied values of n, p, J, and nJ.")
    }
    n <- 197
    p <- 98
    J <- 4
    # Admixture row order: AA (47), EU (50), JP (50), AF (50)
    nJ <- c(47, 50, 50, 50)
  }

  if (missing(nJ)) nJ <- c(rep(n%/%J, J-1), n%/%J + n%%J)
  if (!missing(nJ) & sum(nJ) != n) stop("sum(nJ) must equal n")

  # Generate appropriate scaling factors depending on if partitioning var(y) or not
  scaleXbeta <- sqrt(eta)
  scaleZgamma <- sqrt((1 - eta) * xi)
  scaleEps <- sqrt((1 - eta) * (1 - xi))

  # Gen X
  X <- sim_ps_x(n, nJ, p, structureX, Fst, inbr, standardizeX,
              structureX_other = structureX_other, sampleCols = sampleCols)

  # Gen (scaled) beta and standardized Xbeta (scale by eta)
  j <- 1:p
  s <- rep(1, times = p)
  b <- (j <= p1)
  b <- b * s
  covX <- stats::var(X) * (n - 1)/n
  beta <- b * scaleXbeta / sqrt(drop(crossprod(b, covX) %*% b)) # this scales Xbeta
  Xbeta <- X%*%beta

  # Gen Z
  mlist <- vector("list", J)
  for (jj in 1:J){
    mlist[[jj]] <- matrix(rep(1, nJ[jj]), ncol = 1)
  }
  Z <- as.matrix(Matrix::bdiag(mlist))

  # Gen (scaled) gamma and standardized Zgamma (scale by eta and xi)
  gamma_unscaled <- sim_environ_eff(structureGamma, J)
  env0 <- Z %*% gamma_unscaled
  env <- env0 - mean(env0)
  env_sd <- sqrt(varp(env))
  env <- env/env_sd * scaleZgamma
  ww <- nJ/sum(nJ)
  gamma <- (gamma_unscaled - stats::weighted.mean(gamma_unscaled, ww))/env_sd * scaleZgamma
  # Gen Y
  y <- Xbeta + env + scaleEps * drop(ncvreg::std(as.matrix(stats::rnorm(n))))

  # adding a manual intercept shift
  # probably figure out new SNR later
  if (!missing(B0)) y <- y + B0

  # name things
  w <- 1 + floor(log10(p))
  vlab <- paste0("V", formatC(1:p, format = "d", width = w, flag = "0"))
  colnames(X) <- names(beta) <- vlab

  # return
  list(X = X, y = y, env = env, beta = beta, Z = Z, gamma = gamma, structureX = structureX)

}


