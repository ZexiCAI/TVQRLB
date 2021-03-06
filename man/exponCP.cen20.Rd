% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{exponCP.cen20}
\alias{exponCP.cen20}
\title{Random Length Biased Right Censoring Data with Time Varying Covariates}
\format{A data frame containing 500 rows and 9 columns
\describe{
  \item{S1}{Numerical values giving the first time varying covariate}
  \item{S2}{Numerical values giving the second time varying covariate}
  \item{Z1}{Numerical values giving the first time invariant instrumental variable}
  \item{Z2}{Numerical values giving the second time invariant instrumental variable}
  \item{W1}{Numerical values giving the first change point}
  \item{W2}{Numerical values giving the second change point}
  \item{A}{Numerical values giving the left truncation time}
  \item{Y}{Numerical values giving the observed survival time}
  \item{cens}{Numerical values giving the censoring status. 1 = Failure, 0 = Censored.}
}}
\source{
Randomly generated data
}
\usage{
exponCP.cen20
}
\description{
A dataset containing randomly generated survival data. There are two exponentially distributed change points, and the censoring rate is approximately 20\%.
}
\keyword{datasets}
