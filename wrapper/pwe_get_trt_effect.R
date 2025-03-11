library(dplyr)
library(MCMCpack)

## treatment effect (estimand): difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

#' function to predict survival probability at t years for treated and untreated groups in a piecewise exponential 
#' (PWE) model under priors other than the propensity score-integrated power prior (PSIPP)
#'
#' @param t                 time
#' @param post.samples      posterior samples of a PWE model under various priors (other than PSIPP), with an attribute
#'                          called 'data' which includes the list of variables specified in the data block of the
#'                          Stan program.
#' @param trt.name          name of treatment indicator in the data
#'
get.surv.prob.pwe <- function(
    t,
    post.samples,
    trt.name = "treatment"
){
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  breaks    = stan.data$breaks
  X         = stan.data$X1

  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  lambda    = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz[", 1:J, "]"), drop=F] )
  )
  
  X.trt             = X
  X.trt[, trt.name] = 1 ## assume all subjects are in the treated group
  X.ctl             = X
  X.ctl[, trt.name] = 0 ## assume all subjects are in the control group
  ## Compute linear predictor
  eta.trt           = tcrossprod(beta, X.trt)
  eta.ctl           = tcrossprod(beta, X.ctl)

  ## interval index for which the time t belongs to
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  ## Compute cumulative baseline hazard at each interval
  cumblhaz = apply(lambda, 1, function(x){
    as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
  })
  cumblhaz = matrix(cumblhaz, nrow = J-1)
  cumblhaz = cbind(0, t(cumblhaz))

  ## Compute cumulative hazard
  cumhaz     = lambda[, interval.id] * (t - breaks[interval.id]) + cumblhaz[, interval.id]
  cumhaz.trt = cumhaz * exp(eta.trt)
  cumhaz.ctl = cumhaz * exp(eta.ctl)

  ## predicted t-year survival probability for each subject in data
  S.trt      = exp( -cumhaz.trt )
  S.ctl      = exp( -cumhaz.ctl )

  ## use Bayesian bootstrap to compute predicted t-year survival probability for subjects in the treatment and control arms
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega.trt  = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, nrow(X.trt)))
  omega.ctl  = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, nrow(X.ctl)))
  surv.trt   = rowSums( omega.trt * S.trt)
  surv.ctl   = rowSums( omega.ctl * S.ctl)
  return(
    list(
      surv.trt = surv.trt,
      surv.ctl = surv.ctl
    )
  )
}


#' function to estimate the treatment effect (difference in 2-year survival probability between treated and untreated),
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for a PWE model under priors other than PSIPP
#'
get.surv.diff.2yr.pwe <- function(
    post.samples,
    trt.name = "treatment"
) {
  surv.list = get.surv.prob.pwe (
    t = 2,
    post.samples = post.samples,
    trt.name = trt.name
  )
  surv.diff  = surv.list$surv.trt - surv.list$surv.ctl
  return(surv.diff)
}


#' function to predict survival probability at t years for treated and untreated groups in a PWE model under PSIPP
#'
#' @param t                 time
#' @param post.samples      posterior samples of a PWE model under PSIPP, with an attribute called 'data' which
#'                          includes the list of variables specified in the data block of the Stan program.
#' @param trt.name          name of treatment indicator in the data
#'
get.surv.prob.pwe.psipp <- function(
    t,
    post.samples,
    trt.name = "treatment"
){
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  K         = stan.data$K
  breaks    = stan.data$breaks

  # coefficient of treatment indicator for each stratum
  betaMat   = suppressWarnings(
    as.matrix( post.samples[, paste0(trt.name, '_stratum_', 1:K), drop=F] )
  )
  lambdaMat = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"), drop=F] )
  )
  
  ## interval index for which the time t belongs to
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  ## compute cumulative baseline hazard at t years for each stratum
  cumblhaz = sapply(1:K, function(k){
    lambda = lambdaMat[, paste0("basehaz", "_stratum_", k, "[", 1:J, "]"), drop = F]

    ## Compute cumulative baseline hazard at each interval
    H = apply(lambda, 1, function(x){
      as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
    })
    H = matrix(H, nrow = J-1)
    H = cbind(0, t(H))
    return(
      lambda[, interval.id] * (t - breaks[interval.id]) + H[, interval.id]
    )
  })

  ## predicted t-year survival probability for untreated for each stratum
  S.ctl = lapply(1:K, function(k){
    H_stratum     = cumblhaz[, k]
    S_ctl_stratum = exp( -H_stratum )
    return(S_ctl_stratum)
  }
  )
  S.ctl = do.call(cbind, S.ctl)

  ## predicted t-year survival probability for untreated for each stratum
  S.trt = lapply(1:K, function(k){
    H_stratum     = cumblhaz[, k]
    beta_stratum  = betaMat[, k]
    S_trt_stratum = exp( -H_stratum * exp(beta_stratum) )
    return(S_trt_stratum)
  }
  )
  S.trt = do.call(cbind, S.trt)

  ## take sample mean across strata
  surv.trt   = rowMeans(S.trt)
  surv.ctl   = rowMeans(S.ctl)
  return(
    list(
      surv.trt = surv.trt,
      surv.ctl = surv.ctl
    )
  )
}


#' function to estimate the treatment effect (difference in 2-year survival probability between treated and untreated),
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for a PWE model under PSIPP
#'
get.surv.diff.2yr.pwe.psipp <- function(
    post.samples,
    trt.name = "treatment"
) {
  surv.list = get.surv.prob.pwe.psipp(
    t = 2,
    post.samples = post.samples,
    trt.name = trt.name
  )
  surv.diff  = surv.list$surv.trt - surv.list$surv.ctl
  return(surv.diff)
}
