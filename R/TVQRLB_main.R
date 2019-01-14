#' Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling
#'
#' Fits a quantile regression model to the provided dataset where the covariates are viewed
#' as time-dependent and the sampling is length-biased. The parameters are obtained by
#' minimizing the Euclidean norm of certain estimating equations. For the standard error
#' estimation, standard multiplier bootstrap method is used.
#'
#' @param dataset The survival data.
#' @param betao The initial estimate for the parameter.
#' @param bootstrap_time The bootstrapping time for multiplier bootstrap.
#' @param qtile The quantile level used to conduct the quantile regression.
#'
#' @import survival
#' @import magrittr
#' @export
#' @examples
#' \donttest{
#' TVQRLB(fixedCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
#' TVQRLB(fixedCP.cen40, c(-1, 0.5, 1.5), 100, 0.5)
#' TVQRLB(fixedCP.cen60, c(-1, 0.5, 1.5), 100, 0.5)
#' TVQRLB(exponCP.cen20, c(-1, 0.5, 1.5), 100, 0.5)
#' TVQRLB(exponCP.cen40, c(-1, 0.5, 1.5), 100, 0.5)
#' TVQRLB(exponCP.cen60, c(-1, 0.5, 1.5), 100, 0.5)
#' }
#'
#' @return This function returns a list of lists with each list containing four elements:
#' \itemize{
#'   \item qtile, the quantile level specified to fit the model
#'   \item beta_est, the estimated value of parameter
#'   \item mean_bootstrap, the mean of bootstrap estimates
#'   \item sd_bootstrap, the standard error of bootstrap estimates
#' }
#'
#' @references Cai, Z. and Sit, T. (2019+),
#' "Quantile regression model with time-varying covariates under length-biased sampling,"
#' \emph{Working Paper}.
#'
#' @references Jin, Z., Lin, D., Wei, L., and Ying, Z. (2003),
#' "Rank-based inference for the accelerated failure time model,"
#' \emph{Biometrika}, \bold{90}, 341-353.

TVQRLB <- function(dataset, betao, bootstrap_time, qtile){

  CPs <- "W" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (CPs == 0){
    stop("There are no change points detected. Make sure you name it as \"W1\", \"W2\", etc.")
  }

  covariates <- "S" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (covariates != CPs){
    stop("Number of time-dependent covariates differs from number of change points. Make sure you name it as \"S1\", \"S2\", etc.")
  }

  instruments <- "Z" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (instruments != CPs | instruments != covariates){
    stop("Number of time-invariant instrumental variables differs from number of change points or time-dependent covariates. Make sure you name it as \"Z1\", \"Z2\", etc.")
  }

  survtime <- "\\<A\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (survtime != 1){
    stop("There are no survival time detected. Make sure you name it exactly as \"Y\".")
  }

  truntime <- "\\<A\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (truntime != 1){
    stop("There are no survival time detected. Make sure you name it exactly as \"A\".")
  }

  censtatus <- "\\<cens\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (censtatus != 1){
    stop("There are no censoring status detected. Make sure you name it exactly as \"cens\".")
  }

  if (length(betao) != 1 + CPs){
    stop("The length of initial guess is not correct. Make sure it is the number of change points plus one.")
  }

  if (length(bootstrap_time)>1){
    warning("Multiple bootstrap times are detected. Only the first element is used.")
    bootstrap_time <- bootstrap_time[1]
  }

  sample_size <- bootstrap_size <- length(dataset[,1])

  if (bootstrap_time < sqrt(sample_size)){
    warning("The bootstrap time is inadequate. Change to 1000 instead.")
    bootstrap_time <- 1000
  }

  for (qt in qtile){
    if (qt < 0 | qt > 1){
      stop("Some quantile level specified is out of range. Make sure it is between 0 and 1.")
    }
  }

  res_lst <- vector("list", length = length(qtile))
  for (i in 1:length(qtile)){

    integral_KM_G <- G_surv(dataset$Y-dataset$A, dataset$cens, dataset$Y, rep(1, sample_size))
    beta_est <- optim(fn = U_equation,
                      par = betao,
                      denom = integral_KM_G,
                      dataset = dataset,
                      qtile = qtile[i],
                      weight = rep(1, sample_size),
                      norm = TRUE)$par

    beta_bootstrap <- matrix(nrow = bootstrap_time, ncol = length(betao))
    for (h in 1:bootstrap_time){
      perturb <- rexp(bootstrap_size, 1)
      perturb <- perturb / mean(perturb)

      integral_KM_G_ptb <- G_surv(dataset$Y-dataset$A,
                                  dataset$cens,
                                  dataset$Y,
                                  perturb)

      beta_bootstrap[h, ] <- optim(fn = U_equation,
                                   par = beta_est,
                                   denom = integral_KM_G_ptb,
                                   qtile = qtile[i],
                                   dataset = dataset,
                                   weight = perturb,
                                   norm = TRUE)$par

    }

    print(paste0("FOR QUANTILE LEVEL ", qtile[i]))
    print(paste0("The parameter estimate is: ", beta_est %>% round(4) %>% paste(collapse=" ")))
    print(paste0("The mean of bootstrap estimate is: ", colMeans(beta_bootstrap) %>% round(4) %>% paste(collapse=" ")))
    print(paste0("The standard deviation of bootstrap esimate is: ", cov(beta_bootstrap) %>% diag() %>% sqrt() %>% round(4) %>% paste(collapse=" ")))

    res_lst[[i]] <- list(qtile = qtile[i],
                         beta_est = beta_est,
                         mean_bootstrap = colMeans(beta_bootstrap),
                         sd_bootstrap = cov(beta_bootstrap) %>% diag() %>% sqrt())
  }
  res_lst %>% invisible()
}



#' Quantile Regression Model with Time-Varying Covariates under Length-Biased Sampling using OP Method for Standard Error Estimation
#'
#' Fits a quantile regression model to the provided dataset where the covariates are viewed
#' as time-dependent and the sampling is length-biased. The parameters are obtained by
#' minimizing the Euclidean norm of certain estimating equations. For the standard error
#' estimation, a more computationally efficient algorithm named "Orthogonal Procrustes"
#' method, based on matrix singular value decomposition, is used.
#'
#' @param dataset The survival data.
#' @param betao The initial estimate for the parameter.
#' @param bootstrap_time The bootstrapping time for estimation of the "meat" matrix V.
#' @param qtile The quantile level used to conduct the quantile regression.
#' @param B_size The resampling size for the OP method. Default is \code{10000}.
#'
#' @import survival
#' @import magrittr
#' @importFrom MASS ginv
#' @importFrom pracma sqrtm
#' @export
#' @examples
#' \donttest{
#' TVQRLB_OP(fixedCP.cen20, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' TVQRLB_OP(fixedCP.cen40, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' TVQRLB_OP(fixedCP.cen60, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' TVQRLB_OP(exponCP.cen20, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' TVQRLB_OP(exponCP.cen40, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' TVQRLB_OP(exponCP.cen60, c(-1, 0.5, 1.5), 1000, 0.5, 10000)
#' }
#'
#' @return This function returns a list of lists with each list containing three elements:
#' \itemize{
#'   \item qtile, the quantile level specified to fit the model
#'   \item beta_est, the estimated value of parameter
#'   \item sd_OP, the standard error of parameter estimates
#' }
#'
#' @references Cai, Z. and Sit, T. (2019+),
#' "Quantile regression model with time-varying covariates under length-biased sampling,"
#' \emph{Working Paper}.
#'

TVQRLB_OP <- function(dataset, betao, bootstrap_time, qtile, B_size = 10000){

  CPs <- "W" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (CPs == 0){
    stop("There are no change points detected. Make sure you name it as \"W1\", \"W2\", etc.")
  }

  covariates <- "S" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (covariates != CPs){
    stop("Number of time-dependent covariates differs from number of change points. Make sure you name it as \"S1\", \"S2\", etc.")
  }

  instruments <- "Z" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (instruments != CPs | instruments != covariates){
    stop("Number of time-invariant instrumental variables differs from number of change points or time-dependent covariates. Make sure you name it as \"Z1\", \"Z2\", etc.")
  }

  survtime <- "\\<A\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (survtime != 1){
    stop("There are no survival time detected. Make sure you name it exactly as \"Y\".")
  }

  truntime <- "\\<A\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (truntime != 1){
    stop("There are no survival time detected. Make sure you name it exactly as \"A\".")
  }

  censtatus <- "\\<cens\\>" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()
  if (censtatus != 1){
    stop("There are no censoring status detected. Make sure you name it exactly as \"cens\".")
  }

  if (length(betao) != 1 + CPs){
    stop("The length of initial guess is not correct. Make sure it is the number of change points plus one.")
  }

  if (length(bootstrap_time)>1){
    warning("Multiple bootstrap times are detected. Only the first element is used.")
    bootstrap_time <- bootstrap_time[1]
  }

  sample_size <- bootstrap_size <- length(dataset[,1])
  if (bootstrap_time < sqrt(sample_size)){
    warning("The bootstrap time is inadequate. Change to 1000 instead.")
    bootstrap_time <- 1000
  }

  if (length(B_size)>1){
    warning("Multiple bootstrap sizes are detected. Only the first element is used.")
    B_size <- B_size[1]
  }

  if (B_size < 1000){
    if (B_size < 3){
      warning("The resampling size is inadequate. Change to 10000 instead.")
      B_size <- 10000
    } else{
      warning("The resampling size may not be adequate and may give inacurate results.")
    }
  }

  for (qt in qtile){
    if (qt < 0 | qt > 1){
      stop("Some quantile level specified is out of range. Make sure it is between 0 and 1.")
    }
  }

  res_lst <- vector("list", length = length(qtile))
  for (i in 1:length(qtile)){


    integral_KM_G <- G_surv(dataset$Y-dataset$A, dataset$cens, dataset$Y, rep(1, sample_size))
    beta_est <- optim(fn = U_equation,
                      par = betao,
                      denom = integral_KM_G,
                      dataset = dataset,
                      qtile = qtile[i],
                      weight = rep(1, sample_size),
                      norm = TRUE)$par

    U_bootstrap <- matrix(nrow = bootstrap_time, ncol = length(beta_est))

    for (h in 1:bootstrap_time){
      perturb <- rexp(bootstrap_size, 1)
      perturb <- perturb/mean(perturb)

      bootstrap_integral_KM_G <- G_surv(dataset$Y-dataset$A,
                                        dataset$cens,
                                        dataset$Y,
                                        perturb)

      U_bootstrap[h,] <- sqrt(sample_size) * U_equation(beta_est,
                                                        dataset = dataset,
                                                        denom = bootstrap_integral_KM_G,
                                                        qtile = qtile[i],
                                                        weight = perturb,
                                                        norm = FALSE)
    }
    V_hat <- cov(U_bootstrap)


    A_hat <- matrix(0,length(beta_est),length(beta_est))
    A_hat_inv <- matrix(0,length(beta_est),length(beta_est))


    scl <- 1/1.6
    pass <- FALSE
    while (pass != TRUE){
      scl <- scl * 1.6
      Z_OP <- matrix(scl*(rnorm(length(beta_est)*B_size,0,1)-0.0), nrow = length(beta_est), byrow = TRUE)
      beta_tilde <- matrix(rep(beta_est, B_size),nrow = length(beta_est)) + sample_size^(-0.5) * Z_OP

      ptb <- rexp(sample_size, 1)
      ptb <- ptb/mean(ptb)
      ptb_integral_KM_G <- G_surv(dataset$Y-dataset$A,dataset$cens,dataset$Y,ptb)

      U_OP <- sqrt(sample_size) * apply(beta_tilde,
                                        2,
                                        U_equation,
                                        dataset = dataset,
                                        denom = ptb_integral_KM_G,
                                        qtile = qtile[i],
                                        weight = ptb,
                                        norm = FALSE)

      compst <- sqrt(sample_size) * U_equation(beta_est,
                                               dataset = dataset,
                                               denom = ptb_integral_KM_G,
                                               qtile = qtile[i],
                                               weight = ptb,
                                               norm = FALSE)

      res <- try({ sqrtmU <- pracma::sqrtm(cov(t(U_OP))/scl^2)$B

      W <- MASS::ginv(sqrtmU) %*% (U_OP - matrix(rep(compst,B_size),nrow=length(beta_est)))

      Q <- matrix(0,length(beta_est),length(beta_est))
      R <- matrix(0,length(beta_est),length(beta_est))

      Q <- W %*% t(Z_OP) / B_size

      Q_svd <- svd(Q)
      R <- Q_svd$u %*% diag(c( rep(1, length(beta_est)-1),det(Q_svd$u)*det(Q_svd$v) )) %*% t(Q_svd$v)

      A_hat <- sqrtmU %*% R
      A_hat_inv <- MASS::ginv(A_hat)
      })

      if(inherits(res, "try-error")){
        print(paste0("Warning: Insufficient perturbation, try another set of perturbation."))
        pass <- FALSE
      } else{
        pass <- TRUE
      }
    }

    sd_bootstrap_OP <- rep(NA, length(beta_est))
    sd_bootstrap_OP <- ( (  A_hat_inv %*% V_hat %*% t( A_hat_inv ) ) / sample_size^(1) ) %>% diag() %>% sqrt()

    print(paste0("FOR QUANTILE LEVEL ", qtile[i]))
    print(paste0("The parameter estimate is: ", beta_est %>% round(4) %>% paste(collapse=" ")))
    print(paste0("The standard deviation using OP method is: ", sd_bootstrap_OP %>% round(4) %>% paste(collapse=" ")))

    res_lst[[i]] <- list(qtile = qtile[i],
                         beta_est = beta_est,
                         sd_OP = sd_bootstrap_OP)
  }
  res_lst %>% invisible()
}




#' Form the Estimating Equation
#'
#' Set up an estimating equation.
#'
#' @param beta The parameter estimate.
#' @param dataset The survival data.
#' @param denom The denominator used in the estimating equation.
#' @param qtile The quantile level used to conduct the quantile regression.
#' @param weight The weight of each observation. Default is equal weights.
#' @param norm Whether the norm of the function value should be returned. \code{TRUE}: return the norm; \code{FALSE}: return the function value. Default is \code{TRUE}.
#'
#' @return This function returns the function value (a vector of length n+1) of the estimating equation evaluated at "beta" if norm = FALSE; otherwise returns the norm of it.
#'
#' @references Cai, Z. and Sit, T. (2019+),
#' "Quantile regression model with time-varying covariates under length-biased sampling,"
#' \emph{Working Paper}.
#'
U_equation <- function(beta, dataset, denom, qtile, weight, norm = TRUE){

  CPs <- "W" %>% grep(dataset %>% colnames(), value=TRUE) %>% length()

  dataset$W0 <- rep(0, length(dataset[,1]))
  dataset$S0 <- rep(0, length(dataset[,1]))
  dataset[ ,paste0("W",CPs+1)] <- 100 * dataset[ ,paste0("W",CPs)]

  #int1 <- dataset$Y * (dataset$Y <= dataset$W1)
  #int2 <- (dataset$W1 + (dataset$Y - dataset$W1)*exp( beta[2]*dataset$S1 )) * (dataset$Y > dataset$W1) * (dataset$Y <= dataset$W2)
  #int3 <- (dataset$W1 + (dataset$W2 - dataset$W1)*exp( beta[2]*dataset$S1 ) + (dataset$Y - dataset$W2)*exp( beta[3]*dataset$S2 )) * (dataset$Y > dataset$W2)

  int <- rep(0, length(dataset[,1]))
  for (i in 1:(CPs+1)){
    intInner <- (dataset$Y - dataset$W0) * (i == 1) + (dataset[ ,paste0("W",1)] - dataset[ ,paste0("W",0)]) * (i != 1)
    for (j in seq(from = 2, to = i-1, length.out = max(i-2,0))){
      intInner <- intInner + (dataset[ ,paste0("W",j)] - dataset[ ,paste0("W",j-1)]) * exp(beta[j] * dataset[ ,paste0("S",j-1)])
    }
    intInner <- intInner + (dataset$Y - dataset[ ,paste0("W",i-1)]) * exp(beta[i] * dataset[ ,paste0("S",i-1)])

    int <- int + intInner * (dataset$Y <= dataset[ ,paste0("W",i)]) * (dataset$Y > dataset[ ,paste0("W",i-1)])
  }
  # integral<- (int1+int2+int3)*exp(beta[1])
  integral<- int*exp(beta[1])

  alpha <- length(dataset[,1])*log(length(dataset[,1])-1)
  indicator_smooth <- 1/(1+exp(alpha*(integral-1)))

  dataset$Z0 <- rep(1, length(dataset[,1]))
  U <- NULL
  for (i in 1:(CPs+1)){
    assign(paste0("U", i), weight*dataset$cens / denom * (indicator_smooth-qtile) * dataset[ ,paste0("Z",i-1)])
    U <- methods::cbind2(U, eval(parse(text = paste0("U",i))))
  }
  # U1 <- weight*dataset$cens / denom * (indicator_smooth-qtile) * 1
  # U2 <- weight*dataset$cens / denom * (indicator_smooth-qtile) * dataset$Z1
  # U3 <- weight*dataset$cens / denom * (indicator_smooth-qtile) * dataset$Z2
  # eqn <- colMeans(methods::cbind2(U1, U2, U3))
  eqn <- colMeans(U)

  dataset$Z0 <- NULL
  dataset$S0 <- NULL
  dataset$W0 <- NULL
  dataset[ ,paste0("W",CPs+1)] <- NULL

  if (norm == TRUE){return(norm(eqn, type = "2"))}
  if (norm != TRUE){return(eqn)}
}
