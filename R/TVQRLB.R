#' TVQRLB: Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling
#'
#' This package is an implementation of the estimation and inference procedure of a quantile
#' regression model which is based on an estimating equation, with time-varying covariates
#' under length-biased sampling scheme. It fits a quantile regression model to the provided
#' dataset where the covariates are viewed as time-dependent and the sampling is length-biased.
#' The parameters are obtained by minimizing the Euclidean norm of certain estimating equations.
#' For the standard error estimation, two inference procedures are used: standard multiplier
#' bootstrap, and a more computationally efficient algorithm named "Orthogonal Procrustes" method,
#' based on matrix singular value decomposition, is used.
#'
#' @docType package
#' @name TVQRLB-package
#' @aliases TVQRLB-package
#' @references Cai, Z. and Sit, T. (2018+),
#' "Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling,"
#' \emph{Working Paper}.
"_PACKAGE"
