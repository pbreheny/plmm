#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' NB: this function is simply a wrapper for plmm_prep -> plmm_fit -> plmm_format
#' @param X Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to `get_data()`. 
#' * Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
#' @param fbm Logical: should X be treated as filebacked? Relevant only when X is a string to be passed to `get_data()`. Defaults to NULL, using the default setttings of `get_data()` to determine whether X should be stored in memory.
#' @param std_needed Logical: does the supplied X need to be standardized? Defaults to FALSE, since `process_plink()` standardizes the design matrix by default. 
#' By default, X will be standardized internally.
#' @param y Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector. 
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for Spenncath.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/Spenncath penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/Spenncath penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100. 
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#'
#' @return A list including the estimated coefficients on the original scale, as well as other model fitting details 
#' 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' # using admix data 
#' fit_admix1 <- plmm(X = admix$X, y = admix$y, std_needed = TRUE)
#' s1 <- summary(fit_admix1, idx = 99)
#' print(s1)
#' plot(fit_admix1)
#' 
#' # using admix data and k = 50 
#' fit_admix2 <- plmm(X = admix$X, y = admix$y, k = 50, std_needed = TRUE)
#' s2 <- summary(fit_admix2, idx = 99)
#' print(s2)
#' 
#' # an example with p > n:
#' fit_admix3 <- plmm(X = admix$X[1:50, ], y = admix$y[1:50])
#' 
#' # now use PLINK data files
#' \dontrun{
#' 
#' # file-backed example
#' plmm(X = "~/tmp_files/penncath_lite", # adjust this line to fit current machine
#'  fbm = TRUE, trace = TRUE)
#' 
#' penncath_mid <- process_plink(prefix = "penncath_mid", dataDir = plink_example(path="penncath_mid.fam", parent=T))
#' penncath_clinical <- read.csv(plink_example(path="penncath_clinical.csv"))
#' # for the sake of illustration, I use a simple mean imputation for the outcome 
#' penncath_clinical$hdl_impute <- ifelse(is.na(penncath_clinical$hdl), mean(penncath_clinical$hdl, na.rm = T), penncath_clinical$hdl)
#' 
#' # fit with no 'k' specified
#' fit_plink1 <- plmm(X = penncath_mid$X, y = penncath_clinical$hdl_impute, trace = TRUE)
#' summary(fit_plink1, idx = 5)
#' # Runs in ~219 seconds (3.65 mins) on my 2015 MacBook Pro
#' 
#' # fit with 'k = 5' specified (so using RSpectra::svds())
#' fit_plink2 <- plmm(X = penncath_mid$X, y = penncath_clinical$hdl_impute, k = 5, trace = TRUE)
#' # Runs in ~44 seconds on my 2015 MacBook Pro
#' summary(fit_plink2, idx = 5);summary(fit_plink2, idx = 95)
#' 
#' 
#' # case where X is an FBM
#' lite <- get_data("../temp_files/penncath_lite", fbm = TRUE)
#' clinical <- read.csv("../temp_files/penncath_clinical.csv")
#' hdl <- ifelse(is.na(clinical$hdl), mean(clinical$hdl, na.rm = TRUE), clinical$hdl)
#' fit_fbm <- plmm(X = lite, y = hdl, k = 1200)
#' }


plmm <- function(X,
                 fbm = NULL,
                 std_needed = FALSE,
                 y = NULL,
                 k = NULL, 
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma,
                 alpha = 1,
                 # lambda.min = ifelse(n>p, 0.001, 0.05),
                 lambda.min,
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max.iter = 10000,
                 convex = TRUE,
                 dfmax = NULL,
                 warn = TRUE,
                 penalty.factor = NULL,
                 init = NULL,
                 trace = FALSE) {
  
# check X types -------------------------------------------------
  ## string with a filepath -----------------------------
  if("character" %in% class(X)){
    if(is.null(fbm)){
      dat <- get_data(path = X, trace = trace)
    } else if (fbm){
      dat <- get_data(path = X, returnX = FALSE, trace = trace)
    } else if (!fbm){
      dat <- get_data(path = X, returnX = TRUE, trace = trace)
    }
    
    std_X <- dat$std_X
    if("FBM.code256" %in% class(std_X) | "FBM" %in% class(std_X)){
      fbm_flag <- TRUE
    } else {
        fbm_flag <- FALSE
        }
  }
  ## FBM object ----------------
  if("FBM.code256" %in% class(X) | "FBM" %in% class(X)){stop("To analyze data from a file-backed X matrix, meta-data must also be supplied. 
                                       For the 'X' arg, you need to supply either (1) a character string representing the filepath to the .rds object or (2) a list as returned by get_data(). 
                                       If you don't have an .rds object yet, see process_plink() for preparing your data.")}
  ## list -----------------------
  if("list" %in% class(X)){
    # check for X element 
    if(!("std_X" %in% names(X))){stop("The list supplied for the X argument does not have a design matrix element named 'std_X'. Rename as needed and try again.")}
    if("FBM.code256" %in% class(X$std_X) | "FBM" %in% class(X$std_X)){
      dat <- X
      std_X <- dat$std_X
      # designate the dimensions of the standardized design matrix, with only ns columns
      std_X_n <- std_X$nrow
      std_X_p <- std_X$ncol
      fbm_flag <- TRUE
    } else{
      std_X <- obj$std_X
      std_X_n <- std_X$nrow
      std_X_p <- std_X$ncol
      # TODO: add case to handle X passed as list where X is not an FBM
      fbm_flag <- FALSE
    }
    
  }
  # if FBM flag is not 'on' by now, set it 'off'
  if(!exists('fbm_flag')){fbm_flag <- FALSE}

  ## matrix -------------------------------
  if(!fbm_flag){
    if (!inherits(X, "matrix")) {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(X)=="integer") storage.mode(X) <- "double"
    if (typeof(X)=="character") stop("if X is a matrix, it must be a numeric matrix", call.=FALSE)
    
    # designate the dimensions of the original design matrix 
    n <- nrow(X)
    p <- ncol(X) 
     
    ### handle standardization ---------------------------------------------
    if(std_needed){
      std_X <- ncvreg::std(X)
      std_X_details <- list(center = attr(std_X, 'center'),
                            scale = attr(std_X, 'scale'),
                            ns = attr(std_X, 'nonsingular'))
      
      # designate the dimensions of the standardized design matrix, with only ns columns
    } 
    
  }
  # designate dimensions of the 
  std_X_n <- nrow(std_X)
  std_X_p <- ncol(std_X)
  
  #  check y types -------------------------------
  # if y is null, use .fam file 
  if(is.null(y)){
    # default: uses bigSNP naming convention as would be returned in get_data()
    y <- dat$fam$affection 
  }
  
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }  
  
  # set up defaults --------------------------------------------------
  if(is.null(dfmax)){dfmax <- std_X_p + 1}
  
  # default: penalize everything except the intercept, which we will add later
  if(is.null(penalty.factor)){penalty.factor <- rep(1, std_X_p)}
  
  # set default init
  if(is.null(init)){init <- rep(0, std_X_p)}
  
  # coercion for penalty
  penalty <- match.arg(penalty)
  
  # set default gamma
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  # error checking ------------------------------------------------------------
  if(!fbm_flag){
    # error check for matrix X
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
    if (any(is.na(y)) | any(is.na(std_X))) stop("Missing data (NA's) detected.  
                                            \nTake actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty.factor)!=std_X_p) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  } else {
    #  error checking for FBM X 
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
    if (any(is.na(y))) stop("Missing data (NA's) detected in the outcome.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty.factor)!=std_X_p) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  }
  
  # check K types -------------------------------------------------------
  if (!is.null(K)){
    # first, check type/class:
    if (!inherits(K, "matrix") & !is.list(K)) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be either (1) able to be coerced to a matrix or (2) be a list.", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(std_X) <- "double" # change K to X 
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (is.list(K)) {
      if(!('s' %in% names(K) & 'U' %in% names(K))){stop('Components s and U not both found in list supplied for K.')}
    }
    # last check: look at dimensions
    if (is.matrix(K)){
      if (dim(K)[1] != std_X_n || dim(K)[2] != std_X_n) stop("Dimensions of K and X do not match", call.=FALSE)
    } else if (is.list(K)) {
      if (std_X_n != nrow(K$U)) stop("Dimensions of K and X do not match.")
    }
  
  }
  
  # create a list that captures the centering/scaling for std_X; will need this 
  # later, see `untransform()`
  if(fbm_flag){
    std_X_details <- list(
      center = dat$std_X_center,
      scale = dat$std_X_scale,
      ns = dat$ns
    )
  } 
  


  # prep (SVD)-------------------------------------------------
  if(trace){cat("\nInput data passed all checks.")}
  if (fbm_flag){
    the_prep <- plmm_prep(std_X = std_X,
                          std_X_n = std_X_n,
                          std_X_p = std_X_p,
                          n = dat$n,
                          p = dat$p,
                          y = y,
                          K = K,
                          k = k,
                          diag_K = diag_K,
                          fbm_flag = fbm_flag,
                          trace = trace)
  } else {
    the_prep <- plmm_prep(std_X = std_X,
                          std_X_n = std_X_n,
                          std_X_p = std_X_p,
                          n = n,
                          p = p,
                          y = y,
                          K = K,
                          k = k,
                          diag_K = diag_K,
                          fbm_flag = fbm_flag,
                          trace = trace)
  }
  

  # rotate & fit -------------------------------------------------------------
  if(trace){cat("\nBeginning model fitting.\n")}
  
  the_fit <- plmm_fit(prep = the_prep,
                      std_X_details = std_X_details,
                      eta_star = eta_star,
                      penalty.factor = penalty.factor,
                      fbm_flag = fbm_flag,
                      penalty = penalty,
                      gamma = gamma,
                      alpha = alpha,
                      lambda.min = lambda.min,
                      nlambda = nlambda,
                      lambda = lambda,
                      eps = eps,
                      max.iter = max.iter,
                      warn = warn,
                      init = init)

  if(trace){cat("\nBeta values are estimated -- almost done!")}
  
  # format results ---------------------------------------------------
  if(trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}
  if(fbm_flag){
    # TODO: work format through for FBM case
    the_final_product <- plmm_format(fit = the_fit
                                     # convex = convex,
                                     # dfmax = dfmax, 
                                     )
    
  } else {
    the_final_product <- plmm_format(fit = the_fit
                                     # convex = convex,
                                     # dfmax = dfmax, 
                                     )
    
  }


  return(the_final_product)
  
  
}
