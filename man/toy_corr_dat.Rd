% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toy_corr_dat.R
\name{toy_corr_dat}
\alias{toy_corr_dat}
\title{A function to generate correlated data for testing purposes}
\usage{
toy_corr_dat(n = 100, p = 256, s = 4, gamma = 6, beta = 2, B = 20)
}
\arguments{
\item{n}{Number of observations}

\item{p}{Number of features}

\item{s}{Number of significant (truly associated) features}

\item{gamma}{Magnitude of greatest group mean confounding}

\item{beta}{Magnitude of true coefficient values correpsonding to sigificant features}

\item{B}{Number of groups (e.g., batches, families)}
}
\value{
a list with 'toy' data
}
\description{
A function to generate correlated data for testing purposes
}
\keyword{internal}
