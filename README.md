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


## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main branch with all the latest updates 
  
  - `rspectra` is a branch where we are working to make our singular value decomposition approach in `plmm` scale up to GWAS-sized data using [RSpectra](https://github.com/yixuan/RSpectra). 
  
  - `setup-lambda` is a branch where we are resolving a bug in `setup_lambda`. Stay tuned for updates on this. 
  
  - `process-plink` is a branch where we are toying with how to make the package better process [PLINK](https://www.cog-genomics.org/plink/1.9/) files. 
  
  - `blup` is an **archived** branch previously focused on improving the implementation of the Best Linear Unbiased Predictor method 
  
  - `develop/Yujing` is an **archived** branch. We will be deleting this branch soon -- for internal use only. 