
#' Generate a vector of environmental confounding effects
#'
#' This function allows you to simulate an environmental confounding effect. This is used by \code{genDataPS()}, which scales the generated vector in order to comply with the user specified desired variance in the outcome attributable to the environmental effect.
#' @param structureGamma A character argument describing the desired structure.
#' @param J Number of subpopulations.
#' @export

genGammaUnscaled <- function(structureGamma, J){
  UseMethod("genGammaUnscaled")
}

#' @export
genGammaUnscaled.default <- function(structureGamma, J){
  warning("No valid structure specified - defaulting to 'dichotomous_discordant'")
  gamma_unscaled <- 1:J
  return(gamma_unscaled)
}

#' @export
genGammaUnscaled.character <- function(structureGamma = c('dichotomous_discordant',
                                                          'dichotomous_concordant',
                                                          'linear_discordant',
                                                          'linear_concordant',
                                                          'exponential_discordant',
                                                          'exponential_concordant',
                                                          'halfandhalf_discordant',
                                                          'halfandhalf_concordant'), J){
  structureGamma <- match.arg(structureGamma)
  dat <- strsplit(structureGamma, "_")[[1]]
  if (dat[1] == 'dichotomous'){

    if (dat[2] == 'discordant'){
      g <- rep(1, J)
      j <- 1:J
      s <- c(-1, 1)[j%%2 + 1]
      g <- g * s
    } else {
      first <- floor(J/2)
      second <- J - first
      g <- rep(c(1, -1), times = c(first, second))
    }
  } else {
    if (dat[1] == "linear"){
      g <- 1:J
    } else if (dat[1] == "exponential"){
      g <- rep(2, J)^(1:J)
    } else if (dat[1] == 'halfandhalf'){
      first <- floor(J/2)
      second <- J - first
      g <- rep(c(0.5, 1), times = c(first, second))
    }
    g <- rev(g)
    if (dat[2] == "discordant"){
      j <- 1:J
      s <- c(-1, 1)[j%%2 + 1]
      g <- g * s
    }
  }
  return(g)
}


#' @export
genGammaUnscaled.numeric <- function(structureGamma, J){
  if (length(structureGamma) != J){
    stop("Length of numeric argument 'structureGamma' must equal argument 'J'")
  }
  return(structureGamma)
}
