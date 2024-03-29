% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv-plmm.R
\name{cv.plmm}
\alias{cv.plmm}
\title{Cross-validation for plmm}
\usage{
cv.plmm(
  X,
  y,
  k = NULL,
  K = NULL,
  diag_K = NULL,
  eta_star = NULL,
  penalty = c("MCP", "SCAD", "lasso"),
  penalty.factor = rep(1, ncol(X)),
  type = "blup",
  cluster,
  nfolds = 10,
  seed,
  fold = NULL,
  returnY = FALSE,
  returnBiasDetails = FALSE,
  trace = FALSE,
  ...
)
}
\arguments{
\item{X}{Design matrix for model fitting. May include clinical covariates and
other non-SNP data.}

\item{y}{Continuous outcome vector for model fitting.}

\item{k}{An integer specifying the number of singular values to be used in
the approximation of the rotated design matrix. This argument is passed to
\code{RSpectra::svds()}. Defaults to \code{min(n, p) - 1}, where n and p are the dimensions
of the \emph{standardized} design matrix.}

\item{K}{Similarity matrix, in the form of (1) the relatedness matrix estimated
from the data (default), (2) a user-supplied matrix, or (3) a user-supplied
list with components 'd' and 'u' as created by \code{choose_k()}.}

\item{diag_K}{Logical: should K be a diagonal matrix? This would reflect
observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.}

\item{eta_star}{Optional arg. to \code{plmm_prep}. Defaults to NULL.}

\item{penalty}{The penalty to be applied to the model. Either
"MCP" (the default), "SCAD", or "lasso".}

\item{penalty.factor}{Optional arg. to \code{plmm_prep}. Defaults to 1 for
all predictors (except the intercept).}

\item{type}{A character argument indicating what should be returned from
\code{predict.plmm()}. If \code{type == 'lp'} predictions are based on the linear
predictor, \code{$X beta$}. If \code{type == 'blup'} predictions are based on
the \emph{sum} of the linear predictor and the estimated random effect (BLUP).
Defaults to 'blup', as this has shown to be a superior prediction method in
many applications.}

\item{cluster}{\code{cv.plmm()} can be run in parallel across a cluster using the parallel
package. The cluster must be set up in advance using \code{parallel::makeCluster()}.
The cluster must then be passed to \code{cv.plmm()}.}

\item{nfolds}{The number of cross-validation folds. Default is 10.}

\item{seed}{You may set the seed of the random number generator in order to
obtain reproducible results.}

\item{fold}{Which fold each observation belongs to. By default the observations
are randomly assigned.}

\item{returnY}{Should \code{cv.plmm()} return the linear predictors from the
cross-validation folds? Default is FALSE; if TRUE, this will return a matrix
in which the element for row \code{i}, column \code{j} is the fitted value for observation \code{i}
from the fold in which observation \code{i} was excluded from the fit, at the \code{j}th
value of \code{lambda}.}

\item{returnBiasDetails}{Logical: should the cross-validation bias
(numeric value) and loss (n x p matrix) be returned? Defaults to FALSE.}

\item{trace}{If set to TRUE, inform the user of progress by announcing the
beginning of each CV fold. Default is FALSE.}

\item{...}{Additional arguments to \code{plmm_fit}}
}
\value{
a list with 11 items:
\itemize{
\item type: the type of prediction used ('lp' or 'blup')
\item cve: numeric vector with the cross validation error (CVE) at each value of \code{lambda}
\item cvse: numeric vector with the estimated standard error associated with each value of for \code{cve}
\item fold: numeric \code{n} length vector of integers indicating the fold to which each observation was assigned
\item lambda: numeric vector of \code{lambda} values
\item fit: the overall fit of the object, including all predictors; this is a
list as returned by \code{plmm()}
\item min: The index corresponding to the value of \code{lambda} that minimizes \code{cve}
\item lambda.min: The \code{lambda} value at which \code{cve} is minmized
\item min1se: The index corresponding to the value of \code{lambda} within
standard error of that which minimizes \code{cve}
\item lambda1se: largest value of lambda such that error is within 1 standard error of the minimum.
\item null.dev: numeric value representing the deviance for the
intercept-only model. If you have supplied your own \code{lambda} sequence,
this quantity may not be meaningful.
}
}
\description{
Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized
linear mixed models over a grid of values for the regularization parameter \code{lambda}.
}
\examples{
cv_fit <- cv.plmm(X = admix$X, y = admix$y, seed = 321)
\dontrun{
cv_s <- summary.cv.plmm(cv_fit, lambda = "1se")
print(cv_s)
plot(cv_fit)
}




}
