<!-- badges: start -->
[![GitHub version](https://img.shields.io/static/v1?label=GitHub&message=2.2.1&color=blue&logo=github)](https://github.com/pbreheny/plmm)
[![R-CMD-check](https://github.com/pbreheny/plmm/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/plmm/actions)
<!-- badges: end -->

## Welcome 

The `plmmr` package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects. Documentation for this package is in progress. 


## Installation 

To install the latest version of the package: 

```r
devtools::install_github("pbreheny/plmmr")
```

For a description of the motivation of the functions in this package (along with examples) refer to the second module of [this GWAS data tutorial](https://pbreheny.github.io/adv-gwas-tutorial/index.html)

## Latest changes 

The newest features of `plmmr` are: 
  - version 2.2.0: parameter $\eta$ is no longer estimated in each fold of cross-validation

  - version 2.1.0: A new function `mfdr()` for inference on model coefficients. 
  
  - version 2.0.3: An `xgboost` method is now available in `process_plink()`. Check out the documentation for details. This option should be regarded as being in 'beta-testing' mode.
  
## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main branch with all the latest updates

  - `prep_for_cran` is the development branch where we are doing some bug fixes/polishing things up in preparation for our first CRAN submission (coming up soon!)
  
  - `estimate_eta` is a development branch where we worked through our alternative approach for estimating $\eta$. We're keeping this around for reference as we are writing. 
  
  - `fbm` is a development branch where I am working to extend our current methods to analyze data from a design matrix stored as a file-backed object (a Filebacked Big Matrix, or FBM). See package [bigstatsr](https://privefl.github.io/bigstatsr/) for more info on these objects. 
  
  - `sign_flip` is a development branch where I am toying with how to handle the issues caused by +/- signs being flipped as part of truncated SVD. 
  
  - `refine_workflow` is an **archived** branch in which I explored changing the workflow to make cross validation more efficient. This change involved moving the rotation step into `plmm_prep()`, instead of having that step as part of model fitting in `plmm_fit()`. I found that this change was not compatible with cross validation, for reasons that I am currently writing up as part of a paper. Once that paper is done, I will delete this branch. 
