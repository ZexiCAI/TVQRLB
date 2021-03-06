% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TVQRLB_main.R
\name{TVQRLB}
\alias{TVQRLB}
\title{Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling}
\usage{
TVQRLB(dataset, betao, bootstrap_time, qtile)
}
\arguments{
\item{dataset}{The survival data.}

\item{betao}{The initial estimate for the parameter.}

\item{bootstrap_time}{The bootstrapping time for multiplier bootstrap.}

\item{qtile}{The quantile level used to conduct the quantile regression.}
}
\value{
This function returns a list of lists with each list containing four elements:
\itemize{
  \item qtile, the quantile level specified to fit the model
  \item beta_est, the estimated value of parameter
  \item mean_bootstrap, the mean of bootstrap estimates
  \item sd_bootstrap, the standard error of bootstrap estimates
}
}
\description{
Fits a quantile regression model to the provided dataset where the covariates are viewed
as time-dependent and the sampling is length-biased. The parameters are obtained by
minimizing the Euclidean norm of certain estimating equations. For the standard error
estimation, standard multiplier bootstrap method is used.
}
\examples{
\donttest{
TVQRLB(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
TVQRLB(fixedCP.cen40, c(-1, 0.5, 1.5), 100, 0.5)
TVQRLB(fixedCP.cen60, c(-1, 0.5, 1.5), 100, 0.5)
TVQRLB(exponCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
TVQRLB(exponCP.cen40, c(-1, 0.5, 1.5), 100, 0.5)
TVQRLB(exponCP.cen60, c(-1, 0.5, 1.5), 100, 0.5)
}

}
\references{
Cai, Z. and Sit, T. (2018+),
"Quantile regression model with time-varying covariates under length-biased sampling,"
\emph{Working Paper}.

Jin, Z., Lin, D., Wei, L., and Ying, Z. (2003),
"Rank-based inference for the accelerated failure time model,"
\emph{Biometrika}, \bold{90}, 341-353.
}
