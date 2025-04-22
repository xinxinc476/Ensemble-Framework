library(dplyr)
library(MCMCpack)

## treatment effect (estimand): difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

#' function to predict survival probability at t years for treated and untreated groups in a mixture cure rate (CurePWE)
#' model under priors other than the propensity score-integrated power prior (PSIPP).
#' Note that the PWE model for the non-cured population does not include the treatment assignment as a covariate.
#'
#' @param t                 time
#' @param post.samples      posterior samples of a CurePWE model under various priors (other than PSIPP), with an attribute
#'                          called 'data' which includes the list of variables specified in the data block of the
#'                          Stan program.
#'
get.surv.prob.curepwe <- function(
    t,
    post.samples
){
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  breaks    = stan.data$breaks
  X         = stan.data$X1
  
  p_cured   =  suppressWarnings(
    as.numeric(unlist( post.samples[, "p_cured", drop=F] ))
  )
  
  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  lambda    = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz[", 1:J, "]"), drop=F] )
  )
  ## compute linear predictor
  eta       = tcrossprod(beta, X)
  
  ## interval index for which the time t belongs to
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }
  
  ## compute cumulative hazard for the non-cured population
  if( J > 1 ){
    ## compute cumulative baseline hazard at each interval
    cumblhaz = apply(lambda, 1, function(x){
      as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
    })
    cumblhaz = matrix(cumblhaz, nrow = J-1)
    cumblhaz = cbind(0, t(cumblhaz))
    
    cumhaz = lambda[, interval.id] * (t - breaks[interval.id]) + cumblhaz[, interval.id]
    
  }else{
    cumhaz = lambda[, interval.id] * (t - breaks[interval.id])
  }
  cumhaz = cumhaz * exp(eta)
  
  S_noncured = exp( -cumhaz )
  ## predicted t-year survival probability for each subject in data
  S = p_cured + (1 - p_cured) * S_noncured
  
  ## use Bayesian bootstrap to compute predicted t-year survival probability for each subject
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, nrow(X)))
  surv  = rowSums( omega * S )
  
  return(surv)
}


#' function to predict survival probability at t years for treated/untreated group in a CurePWE model under PSIPP.
#' Note that the PWE model for the non-cured population does not include any covariate.
#'
#' @param t                 time
#' @param post.samples      posterior samples of a CurePWE model under PSIPP, with an attribute called 'data' which
#'                          includes the list of variables specified in the data block of the Stan program.
#'
get.surv.prob.curepwe.psipp <- function(
    t,
    post.samples
){
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  K         = stan.data$K
  breaks    = stan.data$breaks
  
  p_curedMat =  suppressWarnings(
    as.matrix( post.samples[, paste0("p_cured_stratum_", 1:K), drop=F] )
  )
  
  # stratum-specific hazard
  lambdaMat  = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"), drop=F] )
  )
  
  ## interval index for which the time t belongs to
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }
  
  ## compute cumulative hazard at t years for each stratum
  if( J > 1 ){
    cumhaz = sapply(1:K, function(k){
      lambda = lambdaMat[, paste0("basehaz", "_stratum_", k, "[", 1:J, "]"), drop = F]
      
      ## compute cumulative baseline hazard at each interval
      cumblhaz = apply(lambda, 1, function(x){
        as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
      })
      cumblhaz = matrix(cumblhaz, nrow = J-1)
      cumblhaz = cbind(0, t(cumblhaz))
      return(
        lambda[, interval.id] * (t - breaks[interval.id]) + cumblhaz[, interval.id]
      )
    })
    
  }else{
    cumhaz = lambdaMat * (t - breaks[interval.id])
  }
  
  S_noncured = exp( -cumhaz )
  
  ## predicted t-year survival probability for each stratum
  S   = sapply(1:K, function(k){
    p_cured = p_curedMat[, k]
    return(p_cured + (1 - p_cured) * S_noncured[, k])
  })
  ## take sample mean across strata
  surv = rowMeans(S)
  
  return(surv)
}
