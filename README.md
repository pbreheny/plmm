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

  - The main model-fitting function `plmm` now has an argument `k` to approximate the decomposition of relationship matrix $\mathbf{K}$. Specifying `k` can **make model fitting faster by orders of magnitude**, with the caveat that choosing a 'best' value `k` is not trivial and will depend on the structure and size of the data. The approximate singular value decomposition is implemented with [RSpectra](https://github.com/yixuan/RSpectra); check out `RSpectra::svds()` for details.  

  - Processing data from [PLINK](https://www.cog-genomics.org/plink/1.9/) files is now done with `bigsnpr` functions (see `process_plink`). This is faster than the previous version. 
  
## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main branch with all the latest updates

  - `check_size` is a development branch where I am re-thinking how data are saved/return, so that object size is manageable.

  - `choose_k` is a development branch where I am thinking about how to choose the tuning parameter for truncated singular value decomposition. 

  - `fam_lmm` is a development branch where I am thinking about how the framework `plmm()` uses for analyzing relationships can be used in contexts where penalization is not necessary. Example: family-based data which are not high-dimensional. 
  
  - `setup-lambda` is a development branch where we are resolving a bug in `setup_lambda`. Stay tuned for updates on this. 
  
  - `blup` is an **archived** branch previously focused on improving the implementation of the Best Linear Unbiased Predictor method 
  
  - `develop/Yujing` is an **archived** branch. We will be deleting this branch soon -- for internal use only. 
