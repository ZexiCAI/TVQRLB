
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

-   There will be a GitHub repository for installation soon.

Example
-------

This is a basic example with some simulation data:

``` r
library(survival)
library(TVQRLB)
TVQRLB(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
#> [1] "FOR QUANTILE LEVEL 0.5"
#> [1] "The parameter estimate is: -0.8866 0.4668 1.333"
#> [1] "The mean of bootstrap estimate is: -0.8713 0.4126 1.3792"
#> [1] "The standard deviation of bootstrap esimate is: 0.0662 0.0879 0.0904"
```

Reference
---------

Cai, Z. and Sit, T. (2018+), "Quantile regression model with time-varying covariates under length-biased sampling," *Working Paper*.
