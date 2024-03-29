% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plmm.R
\name{plmm}
\alias{plmm}
\title{Fit a linear mixed model with non-convex regularization}
\usage{
plmm(
  X,
  y,
  k = NULL,
  K = NULL,
  diag_K = NULL,
  eta_star = NULL,
  penalty = c("MCP", "SCAD", "lasso"),
  gamma,
  alpha = 1,
  lambda.min,
  nlambda = 100,
  lambda,
  eps = 1e-04,
  max.iter = 10000,
  convex = TRUE,
  dfmax = ncol(X) + 1,
  warn = TRUE,
  penalty.factor = rep(1, ncol(X)),
  init = rep(0, ncol(X)),
  trace = FALSE
)
}
\arguments{
\item{X}{Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to \code{get_data()}.
\itemize{
\item Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
}}

\item{y}{Continuous outcome vector. Logistic regression modeling is still in development.}

\item{k}{An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to \code{RSpectra::svds()}. Defaults to \code{min(n, p) - 1}, where n and p are the dimensions of the \emph{standardized} design matrix.}

\item{K}{Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().}

\item{diag_K}{Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.
Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.}

\item{eta_star}{Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.}

\item{penalty}{The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".}

\item{gamma}{The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for Spenncath.}

\item{alpha}{Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/Spenncath penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/Spenncath penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.}

\item{lambda.min}{The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.}

\item{nlambda}{Length of the sequence of lambda. Default is 100.}

\item{lambda}{A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.}

\item{eps}{Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.}

\item{max.iter}{Maximum number of iterations (total across entire path). Default is 10000.}

\item{convex}{Calculate index for which objective function ceases to be locally convex? Default is TRUE.}

\item{dfmax}{Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.}

\item{warn}{Return warning messages for failures to converge and model saturation? Default is TRUE.}

\item{penalty.factor}{A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.}

\item{init}{Initial values for coefficients. Default is 0 for all columns of X.}

\item{trace}{If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.}
}
\value{
A list which includes:
\itemize{
\item \code{beta_vals}: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of \code{lambda}
\item \code{rotated_scale_beta_vals}: the matrix of estimated coefficients on the ~rotated~ scale. This is the scale on which the model was fit.
\item \code{lambda}: a numeric vector of the lasso tuning parameter values used in model fitting.
\item \code{eta}: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
\item \code{s}: a vectof of the eigenvalues of relatedness matrix \code{K}; see \code{relatedness_mat()} for details.
\item \code{U}: a matrix of the eigenvalues of relatedness matrix \code{K}
\item \code{rot_y}: the vector of outcome values on the rotated scale. This is the scale on which the model was fit.
\item \code{linear.predictors}: the matrix resulting from the product of \code{stdrot_X} and the estimated coefficients on the ~rotated~ scale.
\item \code{penalty}: character string indicating the penalty with which the model was fit (e.g., 'MCP')
\item \code{gamma}: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
\item \code{alpha}: numeric value indicating the elastic net tuning parameter.
\item \code{convex.min}: NULL (This is an option we will add in the future!)
\item \code{loss}: vector with the numeric values of the loss at each value of \code{lambda} (calculated on the ~rotated~ scale)
\item \code{penalty.factor}: vector of indicators corresponding to each predictor, where 1 = predictor was penalized.
\item \code{ns_idx}: vector with the indicies of predictors which were constant features (i.e., had no variation).
\item \code{p}: the number of features
\item \code{n}: the number of observations (instances)
\item \code{iter}: numeric vector with the number of iterations needed in model fitting for each value of \code{lambda}
\item \code{converged}: vector of logical values indicating whether the model fitting converged at each value of \code{lambda}
}
}
\description{
This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
NB: this function is simply a wrapper for plmm_prep -> plmm_fit -> plmm_format
}
\examples{
# using admix data 
fit_admix1 <- plmm(X = admix$X, y = admix$y)
s1 <- summary(fit_admix1, idx = 99)
print(s1)
plot(fit_admix1)

# an example with p > n:
fit_admix2 <- plmm(X = admix$X[1:50, ], y = admix$y[1:50])
s2 <- summary(fit_admix2, idx = 99)
print(s2)
plot(fit_admix2) # notice: the default penalty is MCP

}
