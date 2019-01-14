
TVQRLB: Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling
===========================================================================================

The **TVQRLB** package implements the estimation and inference procedure for the quantile regression model which is based on an estimating equation, with time-varying covariates under length-biased sampling scheme.

Installation
------------

You can install the released version of TVQRLB from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("TVQRLB")
library(TVQRLB)
```

You can install the released version of TVQRLB from GitHub with:

``` r
devtools::install_github("ZexiCAI/TVQRLB")
library(TVQRLB)
```

Example
-------

This is a basic example with some simulation data:

``` r
library(survival)
library(TVQRLB)
TVQRLB(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
#> [1] "FOR QUANTILE LEVEL 0.5"
#> [1] "The parameter estimate is: -0.8866 0.4668 1.333"
#> [1] "The mean of bootstrap estimate is: -0.8816 0.4246 1.3742"
#> [1] "The standard deviation of bootstrap esimate is: 0.0792 0.0849 0.0937"
TVQRLB_OP(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5, 10000)
#> [1] "FOR QUANTILE LEVEL 0.5"
#> [1] "The parameter estimate is: -0.8866 0.4668 1.333"
#> [1] "The standard deviation using OP method is: 0.114 0.1197 0.123"
```

Reference
---------

Cai, Z. and Sit, T. (2019+), "Quantile regression model with time-varying covariates under length-biased sampling," *Working Paper*.
