
TVQRLB: Censored Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling
====================================================================================================

The **TVQRLB** package implements the estimation and inference procedure for the censored quantile regression model which is based on an estimating equation, with time-varying covariates under length-biased sampling scheme.

Installation
------------

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
#> [1] "The mean of bootstrap estimate is: -0.8895 0.4237 1.3821"
#> [1] "The standard deviation of bootstrap esimate is: 0.0759 0.0859 0.0857"
TVQRLB_OP(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5, 10000)
#> [1] "FOR QUANTILE LEVEL 0.5"
#> [1] "The parameter estimate is: -0.8866 0.4668 1.333"
#> [1] "The standard deviation using OP method is: 0.1124 0.0594 0.0993"
```

Reference
---------

Cai, Z. and Sit, T. (2020+), "Censored Quantile regression model with time-varying covariates under length-biased sampling," Under Revision.
