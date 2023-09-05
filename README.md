<!-- badges: start -->
[![GitHub version](https://img.shields.io/static/v1?label=GitHub&message=1.0.1&color=blue&logo=github)](https://github.com/areisett/penalizedLMM)
[![R-CMD-check](https://github.com/areisett/penalizedLMM/workflows/R-CMD-check/badge.svg)](https://github.com/areisett/penalizedLMM/actions)
<!-- badges: end -->

## Welcome 

The `penalizedLMM` package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects. Documentation for this package is in progress. 


## Installation 

To install the latest version of the package: 

```r
devtools::install_github("areisett/penalizedLMM")
```

For a description of the motivation of the functions in this package (along with examples) refer to the second module of [this GWAS data tutorial](https://pbreheny.github.io/adv-gwas-tutorial/index.html)

## Latest changes 

The newest features of `penalizedLMM` are: 

  - version 2.1.0: A new function `mfdr()` for inference on model coefficients. 

  - version 2.0.3: An `xgboost` method is now available in `process_plink()`. Check out the documentation for details. This option should be regarded as being in 'beta-testing' mode.
  
## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main branch with all the latest updates

  - `fbm` is a development branch where I am working to extend our current methods to analyze data from a design matrix stored as a file-backed object (a Filebacked Big Matrix, or FBM). See package [bigstatsr](https://privefl.github.io/bigstatsr/) for more info on these objects. 

  - `mfdr` is a development branch where I am working to create an `mfdr()` method for `plmm` objects. See `ncvreg::mfdr()` for inspiration. 

  - `fam_lmm` is a development branch where I am thinking about how the framework `plmm()` uses for analyzing relationships can be used in contexts where penalization is not necessary. Example: family-based data which are not high-dimensional. 

  - `setup-lambda` is a development branch where we are resolving a bug in `setup_lambda`. Stay tuned for updates on this. 
  
  - `blup` is an **archived** branch previously focused on improving the implementation of the Best Linear Unbiased Predictor method 
  
  - `develop/Yujing` is an **archived** branch. We will be deleting this branch soon -- for internal use only. 
