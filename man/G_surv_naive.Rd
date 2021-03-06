% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/G_surv.R
\name{G_surv_naive}
\alias{G_surv_naive}
\title{Survival Function of Residual Censoring Variable}
\usage{
G_surv_naive(time, cens, time_Y, weight)
}
\arguments{
\item{time}{The follow up time for right censored data.}

\item{cens}{The censoring status, 0 for censoring, 1 for failure.}

\item{time_Y}{The upper bound of the integral.}

\item{weight}{The weight for each observation.}
}
\value{
This function returns the function value of the residual censoring variable at a user-specified upper bound.
}
\description{
Evaluate the survival function of residual censoring variable at a user-specified upper bound.
}
