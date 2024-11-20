#' wrapper functions for fitting an AFT model under PSIPP

library(cmdstanr)
library(posterior)
library(psrwe)
library(hdbayes)

#' Estimate propensity scores (PS), assign subjects into different strata based on PS, and obtain stratum-specific power
#' prior parameter values (a0s) using psrwe package
#'
#' @param ps_fml_covs     right side of the propensity score (PS) formula involving the covariates, e.g., ~ sex + age.
#' @param v_arm           column name corresponding to arm (treatment v.s. control) assignment.
#' @param borrow_ctl_only whether to borrow historical information on control arm only. If FALSE, then historical information on
#'                        both the treatment and control arm are borrowed, and the strata are formed based on the PS of all
#'                        subjects in current study. Otherwise, only information on control arm will be borrowed, and the strata
#'                        will be formed based on the PS of subjects in the current study control arm. Defaults to TRUE.
#'
get.strata.data = function(
    data.list,
    ps_fml_covs,
    ps_method       = c("logistic", "randomforest"),
    v_arm           = "treatment",
    ctl_arm_level   = 0,
    borrow_ctl_only = TRUE,
    nstrata         = 5,
    trim_ab         = c("both", "above", "below", "none"),
    total_borrow    = NULL,
    method          = c("distance", "inverse_distance", "n_current", "n_external")
) {
  data             = data.list[[1]]
  ## stack all historical data sets into one historical data set
  histdata         = do.call(rbind, data.list[-1])
  if ( borrow_ctl_only ){
    histdata = histdata[ histdata[, v_arm] == ctl_arm_level, ]
  }
  n              = nrow(data) ## current data sample size
  n0             = nrow(histdata) ## historical data sample size
  data.all       = as.data.frame( rbind(data, histdata) )
  data.all$study = rep(c(1, 0), times = c(n, n0))
  
  ## estimate propensity scores (PS) via psrwe::psrwe_est() function
  data_ps        = psrwe::psrwe_est(
    data = data.all,
    ps_fml = as.formula( paste0("study", ps_fml_covs) ),
    ps_method = ps_method,
    v_grp = "study",
    cur_grp_level = 1,
    v_arm = v_arm,
    ctl_arm_level = ctl_arm_level,
    stra_ctl_only = borrow_ctl_only,
    nstrata = nstrata,
    trim_ab = trim_ab
  )
  
  ## split the number of subjects to be borrowed among strata based on PS via psrwe::psrwe_borrow() function
  if ( is.null(total_borrow) )
    total_borrow = n0
  ps_bor     = psrwe::psrwe_borrow(
    dtaps = data_ps,
    total_borrow = total_borrow,
    method = method
  )
  a0.strata = ps_bor$Borrow$Alpha ## a0 for each stratum
  a0.strata[is.na(a0.strata)] = 0
  
  ## obtain data sets with strata assignments
  data.all.strata = ps_bor$data
  ## remove subjects in historical study who are trimmed based on PS
  data.all.strata = data.all.strata[!is.na(data.all.strata$`_strata_`), ]
  strata          = as.integer( data.all.strata$`_strata_` )
  
  data.list.new   = list(
    curr = data.all.strata[data.all.strata$study == 1, colnames(data)],
    hist = data.all.strata[data.all.strata$study == 0, colnames(histdata)]
  )
  strata.list     = list(
    curr = strata[ data.all.strata$study == 1 ],
    hist = strata[ data.all.strata$study == 0 ]
  )
  res             = list(
    "data.list"        = data.list.new,
    "strata.list"      = strata.list,
    "a0.strata"        = a0.strata,
    "res.psrwe.est"    = data_ps,
    "res.psrwe.borrow" = ps_bor
  )
  return(res)
}


#' function to get Stan data for PSIPP
get.aft.stan.data.psipp = function(
    formula,  ## the RHS typically only includes treatment indicator
    data.list,
    strata.list,
    a0.strata,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE
) {
  ## current data
  data          = data.list[[1]]
  data$stratum  = strata.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  ## re-order data by event indicator and stratum
  data          = data[order(data[, eventind.name], data$stratum), ]
  t             = data[, time.name]
  y             = log(t)
  eventind      = as.integer( data[, eventind.name] )
  X             = stats::model.matrix(formula, data)
  p             = ncol(X)
  
  ## historical data
  histdata         = do.call(rbind, data.list[-1])
  histdata$stratum = strata.list[[2]]
  ## re-order histdata by event indicator and stratum
  histdata         = histdata[order(histdata[, eventind.name], histdata$stratum), ]
  t0               = histdata[, time.name]
  y0               = log(t0)
  eventind0        = as.integer( histdata[, eventind.name] )
  X0               = stats::model.matrix(formula, histdata)
  
  ## get the number of strata
  K = as.integer( max( data$stratum, histdata$stratum ) )
  
  ## get strata assignment for current and historical data
  stratumID.obs  = data$stratum[which(eventind == 1)]
  stratumID.cen  = data$stratum[which(eventind == 0)]
  stratumID0.obs = histdata$stratum[which(eventind0 == 1)]
  stratumID0.cen = histdata$stratum[which(eventind0 == 0)]
  
  ## check a0.strata values
  if ( !( is.vector(a0.strata) & (length(a0.strata) %in% c(1, K)) ) )
    stop("a0.strata must be a scalar or a vector of length ", K)
  a0.strata <- hdbayes:::to.vector(param = a0.strata, len = K)
  if ( any(a0.strata < 0 | a0.strata > 1 ) )
    stop("Each element of a0.strata must be a scalar between 0 and 1")
  
  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = hdbayes:::to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = hdbayes:::to.vector(param = beta.sd, default.value = 10, len = p)
  
  ## Default half-normal prior on scale parameter is N^{+}(0, 10^2)
  if ( !is.null(scale.mean) ){
    if ( !( is.vector(scale.mean) & (length(scale.mean) == 1) ) )
      stop("scale.mean must be a scalar if scale.mean is not NULL")
  }
  scale.mean = hdbayes:::to.vector(param = scale.mean, default.value = 0, len = 1)
  if ( !is.null(scale.sd) ){
    if ( !( is.vector(scale.sd) & (length(scale.sd) == 1) ) )
      stop("scale.sd must be a scalar if scale.sd is not NULL")
  }
  scale.sd = hdbayes:::to.vector(param = scale.sd, default.value = 10, len = 1)
  
  standat = list(
    'dist'            = hdbayes:::dist.to.integer(dist),
    'n'               = length(eventind),
    'n_obs'           = sum(eventind),
    'n_cen'           = sum(1 - eventind),
    'n0_obs'          = sum(eventind0),
    'n0_cen'          = sum(1 - eventind0),
    'p'               = p,
    'y_obs'           = y[which(eventind == 1)],
    'y_cen'           = y[which(eventind == 0)],
    'X_obs'           = X[which(eventind == 1), ],
    'X_cen'           = X[which(eventind == 0), ],
    'y0_obs'          = y0[which(eventind0 == 1)],
    'y0_cen'          = y0[which(eventind0 == 0)],
    'X0_obs'          = X0[which(eventind0 == 1), ],
    'X0_cen'          = X0[which(eventind0 == 0), ],
    'K'               = K,
    'stratumID_obs'   = stratumID.obs,
    'stratumID_cen'   = stratumID.cen,
    'stratumID0_obs'  = stratumID0.obs,
    'stratumID0_cen'  = stratumID0.cen,
    'a0s'             = a0.strata,
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' obtain posterior samples from an AFT model under PSIPP
#'
aft.psipp = function(
    formula,
    data.list,
    strata.list,
    a0.strata,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  # get Stan data for PSIPP
  standat = get.aft.stan.data.psipp(
    formula     = formula,
    data.list   = data.list,
    strata.list = strata.list,
    a0.strata   = a0.strata,
    dist        = dist,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    scale.mean  = scale.mean,
    scale.sd    = scale.sd,
    get.loglik  = get.loglik
  )
  
  aft_psipp     = cmdstanr::cmdstan_model("Stan_aft/aft_psipp.stan")
  
  ## fit model in cmdstanr
  fit = aft_psipp$sample(data = standat,
                         iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                         ...)
  
  ## rename parameters
  p        = standat$p
  X        = standat$X_obs
  K        = standat$K
  oldnames = c( paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"), 
                paste0("scaleVec[", 1:K, "]") )
  newnames = c( paste0( colnames(X), '_stratum_', rep(1:K, each = p) ),
                paste0("scale_stratum_", 1:K))
  
  d = hdbayes:::rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}


#' Estimate the logarithm of the normalizing constant for propensity score-integrated 
#' power prior (PSIPP)
#'
aft.psipp.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X0_obs
  K         = stan.data$K
  oldnames = c( paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"), 
                paste0("scaleVec[", 1:K, "]") )
  newnames = c( paste0( colnames(X), '_stratum_', rep(1:K, each = p) ),
                paste0("scale_stratum_", 1:K))
  colnames(d)[colnames(d) %in% newnames] = oldnames
  d = d[, oldnames, drop=F]
  
  ## compute log normalizing constants (lognc) for half-normal prior on scale
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)
  stan.data$is_prior          = is.prior
  
  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta       = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta       = matrix(beta, nrow = p, ncol = K)
    scale      = as.numeric( pars[paste0("scaleVec[", 1:K, "]")] )
    prior_lp   = sum( sapply(1:K, function(k){
      sum(dnorm(beta[, k], mean = data$beta_mean, sd = data$beta_sd, log = T)) +
        dnorm(scale[k], mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc
    })
    )
    
    Eta0_obs       = data$X0_obs %*% beta
    Eta0_cen       = data$X0_cen %*% beta
    stratumID0_obs = data$stratumID0_obs
    stratumID0_cen = data$stratumID0_cen
    y0_obs         = data$y0_obs
    y0_cen         = data$y0_cen
    a0s            = data$a0s
    
    eta0_obs = sapply(1:length(stratumID0_obs), function(i){
      Eta0_obs[i, stratumID0_obs[i]]
    })
    eta0_cen = sapply(1:length(stratumID0_cen), function(i){
      Eta0_cen[i, stratumID0_cen[i]]
    })
    data_lp = sum( a0s[stratumID0_obs] * hdbayes:::aft_model_obs_lpdf(y0_obs, eta0_obs, scale[stratumID0_obs], data$dist) ) +
      sum( a0s[stratumID0_cen] * hdbayes:::aft_model_cen_lpdf(y0_cen, eta0_cen, scale[stratumID0_cen], data$dist) )
 
    if( !data$is_prior ){
      Eta_obs        = data$X_obs %*% beta
      Eta_cen        = data$X_cen %*% beta
      stratumID_obs  = data$stratumID_obs
      stratumID_cen  = data$stratumID_cen
      y_obs          = data$y_obs
      y_cen          = data$y_cen
      
      eta_obs        = sapply(1:length(stratumID_obs), function(i){
        Eta_obs[i, stratumID_obs[i]]
      })
      eta_cen        = sapply(1:length(stratumID_cen), function(i){
        Eta_cen[i, stratumID_cen[i]]
      })
      data_lp        = data_lp + sum( hdbayes:::aft_model_obs_lpdf(y_obs, eta_obs, scale[stratumID_obs], data$dist) ) +
        sum( hdbayes:::aft_model_cen_lpdf(y_cen, eta_cen, scale[stratumID_cen], data$dist) )
    }
    return(data_lp + prior_lp)
  }
  
  lb = c(rep(-Inf, p*K), rep(0, K))
  ub = rep(Inf, (p+1)*K)
  names(ub) = colnames(d)
  names(lb) = names(ub)
  
  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = d,
        'log_posterior' = log_density,
        'data' = stan.data,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )
  
  ## Return a list of lognc and output from bridgesampling::bridge_sampler
  res = list(
    'lognc'        = bs$logml,
    'bs'           = bs
  )
  return(res)
}


#' Log marginal likelihood of an accelerated failure time (AFT) model under the propensity 
#' score-integrated power prior (PSIPP) (compute lognc of PSIPP via aft.psipp.lognc())

aft.logml.psipp = function(
    post.samples,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')
  
  ## computing log normalizing constant for PSIPP using all data sets
  res.all = aft.psipp.lognc(
    post.samples   = post.samples,
    is.prior       = FALSE,
    bridge.args    = bridge.args
  )
  
  ## sample from psipp
  aft_psipp_prior = cmdstanr::cmdstan_model("Stan_aft/aft_psipp_prior.stan")
  
  fit = aft_psipp_prior$sample(data = stan.data,
                            iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                            ...)
  summ = posterior::summarise_draws(fit)
  
  hist.post.samples = fit$draws(format = 'draws_df')
  attr(x = hist.post.samples, which = 'data') = stan.data
  
  ## compute log normalizing constant for PP using historical data set
  res.hist = aft.psipp.lognc(
    post.samples   = hist.post.samples,
    is.prior       = TRUE,
    bridge.args    = bridge.args
  )
  
  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## the minimum estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat
  res = list(
    'model'        = "PSIPP",
    'logml'        = res.all$lognc - res.hist$lognc,
    'bs'           = res.all$bs,
    'bs.hist'      = res.hist$bs,
    'min_ess_bulk' = min(summ[, 'ess_bulk']),
    'max_Rhat'     = max(summ[, 'rhat'])
  )
  
  if ( res[['min_ess_bulk']] < 1000 )
    warning(
      paste0(
        'The minimum bulk effective sample size of the MCMC sampling is ',
        round(res[['min_ess_bulk']], 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res[['max_Rhat']] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        round(res[['max_Rhat']], 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}


#' Log marginal likelihood of an accelerated failure time (AFT) model under the propensity 
#' score-integrated power prior (PSIPP) (compute lognc of PSIPP via computing sum pf lognc
#' of PP within each stratum)
#'
#aft.logml.psipp = function(
#    post.samples,
#    bridge.args   = NULL,
#    iter_warmup   = 1000,
#    iter_sampling = 1000,
#    chains        = 4,
#    ...
#) {
#  stan.data = attr(post.samples, 'data')
  
#  d        = as.matrix(post.samples)
  ## rename parameters
#  p         = stan.data$p
#  X         = stan.data$X0_obs
#  K         = stan.data$K
#  oldnames = c( paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"), 
#                paste0("scaleVec[", 1:K, "]") )
#  newnames = c( paste0( colnames(X), '_stratum_', rep(1:K, each = p) ),
#                paste0("scale_stratum_", 1:K))
#  colnames(d)[colnames(d) %in% newnames] = oldnames
#  d = d[, oldnames, drop=F]
  
  ## compute log normalizing constants (lognc) for half-normal prior on scale
#  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)
  
  ## compute lognc for PSIPP (equivalent to sum of lognc under PP from each stratum)
#  start.idx0.obs = stan.data$start_idx0_obs
#  end.idx0.obs   = stan.data$end_idx0_obs
#  start.idx0.cen = stan.data$start_idx0_cen
#  end.idx0.cen   = stan.data$end_idx0_cen
#  n0.obs.stratum = end.idx0.obs - start.idx0.obs + 1 ## sample size for each stratum in historical data (uncensored)
#  n0.cen.stratum = end.idx0.cen - start.idx0.cen + 1 ## sample size for each stratum in historical data (censored)
  
#  X0.obs         = stan.data$X0_obs
#  X0.cen         = stan.data$X0_cen
#  y0.obs         = stan.data$y0_obs
#  y0.cen         = stan.data$y0_cen
  
#  aft_pp_prior = cmdstanr::cmdstan_model("Stan_aft/aft_pp_prior.stan")
  
#  res.lognc.pp = lapply(1:K, function(k){
    ## get Stan data for PP using historical data from each stratum
#    hist.stan.data = list(
#      'dist'            = stan.data$dist,
#      'n0_obs'          = n0.obs.stratum[k],
#      'n0_cen'          = n0.cen.stratum[k],
#      'p'               = stan.data$p,
#      'y0_obs'          = y0.obs[ start.idx0.obs[k]:end.idx0.obs[k] ],
#      'y0_cen'          = y0.cen[ start.idx0.cen[k]:end.idx0.cen[k] ],
#      'X0_obs'          = X0.obs[start.idx0.obs[k]:end.idx0.obs[k], ],
#      'X0_cen'          = X0.cen[start.idx0.cen[k]:end.idx0.cen[k], ],
#      'a0'              = stan.data$a0s[k],
#      'beta_mean'       = stan.data$beta_mean,
#      'beta_sd'         = stan.data$beta_sd,
#      'scale_mean'      = stan.data$scale_mean,
#      'scale_sd'        = stan.data$scale_sd
#    )
    ## fit PP using historical data sets
#    fit    = aft_pp_prior$sample(data = hist.stan.data,
#                           iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
#                           ...)
    
#    hist.post.samples = fit$draws(format = 'draws_df')
#    attr(x = hist.post.samples, which = 'data') = hist.stan.data
    
    ## compute log normalizing constant for PP using historical data sets
#    res.hist = hdbayes:::aft.pp.lognc(
#      post.samples   = hist.post.samples,
#      is.prior       = TRUE,
#      bridge.args    = bridge.args
#    )
#    return(res.hist)
#  })
#  lognc_psipp = sum( unlist( lapply(res.lognc.pp, function(l){l[['lognc']]}) ) )
  
  ## log of the unnormalized posterior density function
#  log_density = function(pars, data){
#    p          = data$p
#    K          = data$K
#    beta       = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
#    beta       = matrix(beta, nrow = p, ncol = K)
#    scale      = as.numeric( pars[paste0("scaleVec[", 1:K, "]")] )
#    prior_lp   = sum( sapply(1:K, function(k){
#      sum(dnorm(beta[, k], mean = data$beta_mean, sd = data$beta_sd, log = T)) +
#        dnorm(scale[k], mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc
#    })
#    )
    
#    Eta_obs        = data$X_obs %*% beta
#    Eta_cen        = data$X_cen %*% beta
#    Eta0_obs       = data$X0_obs %*% beta
#    Eta0_cen       = data$X0_cen %*% beta
#    start_idx_obs  = data$start_idx_obs
#    end_idx_obs    = data$end_idx_obs
#    start_idx_cen  = data$start_idx_cen
#    end_idx_cen    = data$end_idx_cen
#    start_idx0_obs = data$start_idx0_obs
#    end_idx0_obs   = data$end_idx0_obs
#    start_idx0_cen = data$start_idx0_cen
#    end_idx0_cen   = data$end_idx0_cen
#    y_obs          = data$y_obs
#    y_cen          = data$y_cen
#    y0_obs         = data$y0_obs
#    y0_cen         = data$y0_cen
#    a0s            = data$a0s
    
#    data_lp    = sum( sapply(1:K, function(k){
#      hdbayes:::aft_model_lp(y_obs[ start_idx_obs[k]:end_idx_obs[k] ], y_cen[ start_idx_cen[k]:end_idx_cen[k] ], 
#                             Eta_obs[start_idx_obs[k]:end_idx_obs[k], k], Eta_cen[start_idx_cen[k]:end_idx_cen[k], k], scale[k], data$dist) + 
#        a0s[k] * hdbayes:::aft_model_lp(y0_obs[ start_idx0_obs[k]:end_idx0_obs[k] ], y0_cen[ start_idx0_cen[k]:end_idx0_cen[k] ], 
#                                        Eta0_obs[start_idx0_obs[k]:end_idx0_obs[k], k], Eta0_cen[start_idx0_cen[k]:end_idx0_cen[k], k], scale[k], data$dist)
#    }) )
#    return(data_lp + prior_lp)
#  }
  
#  lb = c(rep(-Inf, p*K), rep(0, K))
#  ub = rep(Inf, (p+1)*K)
#  names(ub) = colnames(d)
#  names(lb) = names(ub)
  
#  bs = do.call(
#    what = bridgesampling::bridge_sampler,
#    args = append(
#      list(
#        "samples" = d,
#        'log_posterior' = log_density,
#        'data' = stan.data,
#        'lb' = lb,
#        'ub' = ub),
#      bridge.args
#    )
#  )
  
  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## and outputs from computing the log normalizing constant for PSIPP
#  res = list(
#    'model'        = "PSIPP",
#    'logml'        = bs$logml - lognc_psipp,
#    'bs'           = bs,
#    'res.lognc.pp' = res.lognc.pp
#  )
#  return(res)
#}
