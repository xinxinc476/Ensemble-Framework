# for analysis using AFT model
library(dplyr)
library(MCMCpack)
#library(collapse)

## treatment effect (estimand): difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

#' function to predict survival probability at t years for treated and untreated groups for priors other than PSIPP
#'
#' @param t                 time
#' @param post.samples      posterior samples of an AFT model under various priors (other than PSIPP), with an attribute
#'                          called 'data' which includes the list of variables specified in the data block of the
#'                          Stan program.
#' @param trt.name          name of treatment indicator in the data
#'
get.surv.prob <- function(
    t,
    post.samples,
    trt.name = "treatment"
){
  stan.data = attr(post.samples, 'data')
  dist      = stan.data$dist
  X.obs     = stan.data$X_obs
  X.cen     = stan.data$X_cen
  X         = rbind(X.obs, X.cen)
  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  scale     = suppressWarnings(
    as.numeric( unlist(post.samples[, "scale"]) )
  )
  
  eta.mtx    = tcrossprod(beta, X)
  ## predicted t-year survival probability for each subject in data
  S          = apply(eta.mtx, 2, function(c){
    if ( dist == 1 ){
      return(
        stats::pnorm(log(t), mean = c, sd = scale, lower.tail = F)
      )
    }else if( dist == 2 ){
      return(
        stats::plogis(log(t), location = c, scale = scale, lower.tail = F)
      )
    }else if( dist == 3 ){
      z = (-log(t) + c) / scale
      return( exp(-exp(-z)) )
    }
  })
  
  id.trt     = which( X[, trt.name] == 1 )
  id.ctl     = which( X[, trt.name] == 0 )
  
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
get.surv.diff.2yr <- function(
    post.samples,
    trt.name = "treatment"
) {
  surv.list = get.surv.prob(
    t = 2,
    post.samples = post.samples,
    trt.name = trt.name
  )
  surv.diff  = surv.list$surv.trt - surv.list$surv.ctl
  return(surv.diff)
}


#' function to predict survival probability at t years for treated and untreated groups for PSIPP
#'
#' @param t                 time
#' @param post.samples      posterior samples of an AFT model under PSIPP, with an attribute called 'data' which 
#'                          includes the list of variables specified in the data block of the Stan program.
#' @param trt.name          name of treatment indicator in the data                          
#'
get.surv.prob.psipp <- function(
    t,
    post.samples,
    trt.name = "treatment"
){
  stan.data = attr(post.samples, 'data')
  K         = stan.data$K # number of strata
  dist      = stan.data$dist
  
  # intercept for each stratum
  beta0     = suppressWarnings(
    as.matrix( post.samples[, paste0('(Intercept)_stratum_', 1:K), drop=F] )
  )
  # coefficient of treatment indicator for each stratum
  betat     = suppressWarnings(
    as.matrix( post.samples[, paste0(trt.name, '_stratum_', 1:K), drop=F] )
  )
  scale     = suppressWarnings(
    as.matrix( post.samples[, paste0("scale_stratum_", 1:K), drop=F] )
  )
  
  # predicted t-year survival probability for untreated for each stratum
  S.ctl = lapply(1:K, function(k){
    eta_stratum   = as.numeric( beta0[, k] )
    scale_stratum = as.numeric( scale[, k] )
    
    if ( dist == 1 ){
      return(
        stats::pnorm(log(t), mean = eta_stratum, sd = scale_stratum, lower.tail = F)
      )
    }else if( dist == 2 ){
      return(
        stats::plogis(log(t), location = eta_stratum, scale = scale_stratum, lower.tail = F)
      )
    }else if( dist == 3 ){
      z = (-log(t) + eta_stratum) / scale_stratum
      return( exp(-exp(-z)) )
    }
  })
  S.ctl = do.call(cbind, S.ctl)
  
  # predicted t-year survival probability for treated for each stratum
  S.trt = lapply(1:K, function(k){
    eta_stratum   = as.numeric( beta0[, k] + betat[, k] )
    scale_stratum = as.numeric( scale[, k] )
    
    if ( dist == 1 ){
      return(
        stats::pnorm(log(t), mean = eta_stratum, sd = scale_stratum, lower.tail = F)
      )
    }else if( dist == 2 ){
      return(
        stats::plogis(log(t), location = eta_stratum, scale = scale_stratum, lower.tail = F)
      )
    }else if( dist == 3 ){
      z = (-log(t) + eta_stratum) / scale_stratum
      return( exp(-exp(-z)) )
    }
  })
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
    trt.name = "treatment"
) {
  surv.list = get.surv.prob.psipp(
    t = 2,
    post.samples = post.samples,
    trt.name = trt.name
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
