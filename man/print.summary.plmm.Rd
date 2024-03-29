% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-summary-plmm.R
\name{print.summary.plmm}
\alias{print.summary.plmm}
\title{A function to print the summary of a \code{plmm} model}
\usage{
\method{print}{summary.plmm}(x, ...)
}
\arguments{
\item{x}{A \code{summary.plmm} object}

\item{...}{Not used}
}
\description{
A function to print the summary of a \code{plmm} model
}
\examples{
fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
fit2 <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X), penalty = "SCAD")
s1 <- summary(fit, idx = 97)
s2 <- summary(fit, lambda = fit$lambda[97])
s3 <- summary(fit2, idx = 25)
print(s1)
print(s2)
print(s3)
}
