library(tidyverse)
library(MCMCpack)

## treatment effect defined as log hazard ratio at 2 years for treated v.s. untreated

#' function to compute treatment effect for priors other than PSIPP
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
  log_hazard = beta[, 1:J]
  log_hazard = log_hazard[, 1] + cbind(0, log_hazard[, 2:J])

  ## interval index for which the value 2 belongs to
  interval.id = findInterval(2, c(0, breaks, Inf), left.open = TRUE)
  log_hazard  = log_hazard[, interval.id] ## log baseline hazard at 2 years

  beta        = beta[, -(1:J)]
  X.data      = as.matrix( data[, colnames(beta)] )
  eta.mtx     = tcrossprod(beta, X.data)

  ## predicted 2-year log hazard for each subject in data
  log_hazard_2yr = log_hazard + eta.mtx

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

## for PSIPP
get.log.hazard.ratio.2yr.psipp <- function(
    post.samples,
    res.strata,
    breaks
) {
  J          = length(breaks) + 1
  stan.data  = attr(post.samples, 'data')
  p          = stan.data$p
  K          = stan.data$K
  X          = stan.data$X
  beta       = suppressWarnings(
    as.matrix( post.samples[, paste0( "treatment", '_stratum_', 1:K ), drop=F] )
  )
  log_hazard = suppressWarnings(
    as.matrix( post.samples[, paste0( colnames(X)[-ncol(X)], '_stratum_', rep(1:K, each = p-1) ), drop=F] )
  )
  interval_names = colnames(X)[-ncol(X)]

  ## interval index for which the value 2 belongs to
  interval.id        = findInterval(2, c(0, breaks, Inf), left.open = TRUE)
  log_hazard_stratum = sapply(1:K, function(k){
    log_hazard_stratum = log_hazard[, paste0(interval_names, '_stratum_', k)]
    log_hazard_stratum = log_hazard_stratum[, 1] + cbind(0, log_hazard_stratum[, 2:J])
    return( log_hazard_stratum[, interval.id]) ## log baseline hazard at 2 years
  })

  curr.psipp         = res.strata$data.list$curr
  curr.psipp$stratum = as.integer( res.strata$strata.list[[1]] )
  ## re-order curr.psipp by stratum
  curr.psipp         = curr.psipp[order(curr.psipp$stratum), ]
  ## get starting and ending indices for each strata
  num.obs            = as.numeric( table(curr.psipp$stratum) )
  end.idx            = cumsum(num.obs)
  start.idx          = c(1, end.idx[-K] + 1)
  trt.psipp          = curr.psipp[, "treatment", drop = F]

  ## predicted 2-year log hazard for each subject in data
  log_hazard_2yr = lapply(1:K, function(k){
    beta_stratum = beta[, paste0( 'treatment_stratum_', k), drop = F]
    eta_stratum  = tcrossprod(beta_stratum, trt.psipp[ start.idx[k]:end.idx[k], ])
    return(
      eta_stratum + log_hazard_stratum[, k]
    )
  })
  log_hazard_2yr = do.call(cbind, log_hazard_2yr)

  id.trt     = which( curr.psipp$treatment == 1 )
  id.ctl     = which( curr.psipp$treatment == 0 )

  ## use Bayesian bootstrap to compute predicted 2-year survival probability for subjects in the treatment and control arms
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega.trt      = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.trt)))
  omega.ctl      = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, length(id.ctl)))
  log.hazard.trt = rowSums( omega.trt * log_hazard_2yr[, id.trt])
  log.hazard.ctl = rowSums( omega.ctl * log_hazard_2yr[, id.ctl])
  log.hazard.ratio = log.hazard.trt - log.hazard.ctl
  return(log.hazard.ratio)
}
