#' Integral of Survival Function of Residual Censoring Variable
#' 
#' Integrate the survival function of residual censoring variable from 0 to a user-specified upper bound.
#' 
#' @param time The follow up time for right censored data.
#' @param cens The censoring status, 0 for censoring, 1 for failure.
#' @param time_Y The upper bound of the integral.
#' @param weight The weight for each observation.
#' 
#' @import survival stats
#' @export
#' 
#' @examples 
#' G_surv(time = rexp(100,1), cens = rbinom(100, 1, 0.8), time_Y = 1.5, weight = rep(1, 100))
#' G_surv(time = rexp(100,1), cens = rbinom(100, 1, 0.8), time_Y = 1.5, weight = rexp(100, 1))
#' 
#' @return This function returns the integral of the residual censoring variable from 0 to a user-specified upper bound.
#' 

G_surv <- function(time, cens, time_Y, weight){
  
  KM_G_surv <- survfit(Surv(time - 0.001*cens, 1-cens) ~ 1, weights = weight)
  
  delta_KM_G <- diff(c(0, KM_G_surv$time)) * c(1,KM_G_surv$surv[1:(length(KM_G_surv$time)-1)])
  
  intfun <- approxfun(c(0,KM_G_surv$time,1000*max(KM_G_surv$time)), c(0,cumsum(delta_KM_G),max(cumsum(delta_KM_G))))
  
  integral_G <- intfun(time_Y)
  return(integral_G)
}


#' Survival Function of Residual Censoring Variable
#' 
#' Evaluate the survival function of residual censoring variable at a user-specified upper bound.
#' 
#' @param time The follow up time for right censored data.
#' @param cens The censoring status, 0 for censoring, 1 for failure.
#' @param time_Y The upper bound of the integral.
#' @param weight The weight for each observation.
#' 
#' @import survival stats
#' 
#' @return This function returns the function value of the residual censoring variable at a user-specified upper bound.
#' 

G_surv_naive <- function(time, cens, time_Y, weight){
  KM_G_surv <- survfit(Surv(time - 0.001*cens, 1-cens) ~ 1, weights = weight)
  
  for (i in 1:length(KM_G_surv$surv)){
    if (KM_G_surv$surv[i] == 0){
      KM_G_surv$surv[i] <- KM_G_surv$surv[i-1]
    }
  }
  
  stpfun <- stepfun(KM_G_surv$time, c(1,KM_G_surv$surv))
  return(stpfun(time_Y))
}
