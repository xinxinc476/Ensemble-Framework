library(tidyverse)
library(collapse)
library(MCMCpack)

## treatment effect (estimand): difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

#' function to predict survival probability at t years for treated and untreated groups for priors other than PSIPP
#'
#' @param t                 time
#' @param post.samples      posterior samples of a GLM under various priors (other than PSIPP), with an attribute
#'                          called 'data' which includes the list of variables specified in the data block of the
#'                          Stan program.
#' @param data              time-to-event data (not in counting process format)
#' @param breaks            cut points for time intervals
#'
get.surv.prob.glm <- function(
    t,
    post.samples,
    data,
    breaks
){
  J         = length(breaks) + 1
  stan.data = attr(post.samples, 'data')
  X         = stan.data$X
  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  log_hazard = beta[, 1:J]
  log_hazard = log_hazard[, 1] + cbind(0, log_hazard[, 2:J])
  hazard     = exp(log_hazard) ## baseline hazards

  ## interval index for which the time t belongs to
  interval.id = findInterval(t, c(0, breaks, Inf), left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  breaks.new              = (c(breaks, Inf))[1:interval.id]
  breaks.new[interval.id] = t
  interval.length         = breaks.new - c(0, breaks.new[-interval.id])

  H          = hazard[, 1:interval.id, drop = F] %r*% interval.length
  if( interval.id != 1 ){
    H        = t( apply(H, 1, cumsum) ) ## cumulative baseline hazard
  }
  H          = H[, interval.id] ## cumulative baseline hazard at t years
  beta       = beta[, -(1:J)]

  X.data     = as.matrix( data[, colnames(beta)] )
  eta.mtx    = tcrossprod(beta, X.data)
  ## predicted t-year survival probability for each subject in data
  S          = exp( -H * exp(eta.mtx) )

  id.trt     = which( data$treatment == 1 )
  id.ctl     = which( data$treatment == 0 )

  ## use Bayesian bootstrap to compute predicted t-year survival probability for subjects in the treatment and control arms
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega.trt  = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.trt)))
  omega.ctl  = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.ctl)))
  surv.trt   = rowSums( omega.trt * S[, id.trt])
  surv.ctl   = rowSums( omega.ctl * S[, id.ctl])
  return(
    list(
      surv.trt = surv.trt,
      surv.ctl = surv.ctl
    )
  )
}


#' function to estimate the treatment effect (difference in 2-year survival probability between treated and untreated),
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for priors other than PSIPP
get.surv.diff.2yr.glm <- function(
    post.samples,
    data,
    breaks
) {
  surv.list = get.surv.prob.glm(
    t = 2,
    post.samples = post.samples,
    data = data,
    breaks = breaks
  )
  surv.diff  = surv.list$surv.trt - surv.list$surv.ctl
  return(surv.diff)
}


#' function to predict survival probability at t years for treated and untreated groups for PSIPP
#'
#' @param t                 time
#' @param post.samples      output from [glm.stratified.pp()] giving posterior samples of a GLM under the stratified power
#'                          prior (PP), with an attribute called 'data' which includes the list of variables specified in
#'                          the data block of the Stan program.
#' @param breaks            cut points for time intervals
#'
get.surv.prob.psipp <- function(
    t,
    post.samples,
    breaks
){
  J          = length(breaks) + 1
  stan.data  = attr(post.samples, 'data')
  p          = stan.data$p
  K          = stan.data$K
  X          = stan.data$X
  beta       = suppressWarnings(
    as.matrix( post.samples[, paste0( "treatment", '_stratum_', 1:K ), drop=F] )
  )

  log_hazard = suppressWarnings(
    as.matrix( post.samples[, paste0( colnames(X)[-which(colnames(X) == "treatment")], '_stratum_', rep(1:K, each = p-1) ), drop=F] )
  )
  interval_names = colnames(X)[-which(colnames(X) == "treatment")]

  ## interval index for which the time t belongs to
  interval.id = findInterval(t, c(0, breaks, Inf), left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  breaks.new              = (c(breaks, Inf))[1:interval.id]
  breaks.new[interval.id] = t
  interval.length         = breaks.new - c(0, breaks.new[-interval.id])

  ## compute cumulative baseline hazard at t years for each stratum
  H = sapply(1:K, function(k){
    log_hazard_stratum = log_hazard[, paste0(interval_names, '_stratum_', k)]
    log_hazard_stratum = log_hazard_stratum[, 1] + cbind(0, log_hazard_stratum[, 2:J])
    hazard_stratum     = exp(log_hazard_stratum) ## baseline hazards

    H_stratum          = hazard_stratum[, 1:interval.id, drop = F] %r*% interval.length
    if( interval.id != 1 ){
      H_stratum        = t( apply(H_stratum, 1, cumsum) ) ## cumulative baseline hazard
    }
    H_stratum        = H_stratum[, interval.id]
    return( H_stratum) ## cumulative baseline hazard at 2 years
  })

  ## predicted t-year survival probability for untreated for each stratum
  S.ctl = lapply(1:K, function(k){
    H_stratum     = H[, k]
    S_ctl_stratum = exp( -H_stratum )
    return(S_ctl_stratum)
  }
  )
  S.ctl = do.call(cbind, S.ctl)

  ## predicted t-year survival probability for untreated for each stratum
  S.trt = lapply(1:K, function(k){
    H_stratum     = H[, k]
    beta_stratum  = beta[, paste0( 'treatment_stratum_', k)]
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
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for PSIPP
get.surv.diff.2yr.psipp <- function(
    post.samples,
    breaks
) {
  surv.list = get.surv.prob.psipp(
    t = 2,
    post.samples = post.samples,
    breaks = breaks
  )
  surv.diff  = surv.list$surv.trt - surv.list$surv.ctl
  return(surv.diff)
}


#' function to obtain (posterior) samples from the model-averaged priors based on (posterior) samples
#' from individual priors and (posterior) weights
#'
sample.model.avg.prior = function(wts, samples.mtx){
  wts <- as.numeric(wts)
  ## draw n i.i.d. samples (c0) from categorical distribution with probability being `wts`
  c0 <- sample(x = seq_len(ncol(samples.mtx)), size = nrow(samples.mtx), replace = T,
               prob = wts)

  models <- unique(c0)
  res.samples <- lapply(models, function(j){
    nsample = sum(c0 == j)
    return( sample(samples.mtx[, j], size = nsample, replace = T) )
  })
  res.samples <- unlist(res.samples)
  return(res.samples)
}


#' function to compute log hazard ratio at 2 years for treated v.s. untreated for priors other than PSIPP
#' @param data  generated current data set (time-to-event data, not in counting process format)
get.log.hazard.ratio.2yr.glm <- function(
    post.samples,
    data,
    breaks
) {
  J         = length(breaks) + 1
  stan.data = attr(post.samples, 'data')
  X         = stan.data$X
  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  beta        = beta[, -(1:J)]
  X.data      = as.matrix( data[, colnames(beta)] )
  eta.mtx     = tcrossprod(beta, X.data)

  ## predicted 2-year log hazard for each subject in data - log baseline hazard at 2 years
  ## as both treated and untreated have the same value of log baseline hazard at 2 years
  log_hazard_2yr = eta.mtx

  id.trt     = which( data$treatment == 1 )
  id.ctl     = which( data$treatment == 0 )

  ## use Bayesian bootstrap to compute predicted 2-year log hazard for subjects in the treatment and control arms
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega.trt      = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.trt)))
  omega.ctl      = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.ctl)))
  log.hazard.trt = rowSums( omega.trt * log_hazard_2yr[, id.trt])
  log.hazard.ctl = rowSums( omega.ctl * log_hazard_2yr[, id.ctl])
  log.hazard.ratio = log.hazard.trt - log.hazard.ctl
  return(log.hazard.ratio)
}


#' function to compute log hazard ratio at 2 years for treated v.s. untreated for PSIPP
get.log.hazard.ratio.2yr.psipp <- function(
    post.samples,
    breaks
) {
  J          = length(breaks) + 1
  stan.data  = attr(post.samples, 'data')
  p          = stan.data$p
  K          = stan.data$K
  X          = stan.data$X

  ## estimated log hazard ratio for treated v.s. untreated within each stratum for each time interval
  beta       = suppressWarnings(
    as.matrix( post.samples[, paste0( "treatment", '_stratum_', 1:K ), drop=F] )
  )

  ## take sample mean across strata
  log.hazard.ratio = rowMeans(beta)
  return(log.hazard.ratio)
}
