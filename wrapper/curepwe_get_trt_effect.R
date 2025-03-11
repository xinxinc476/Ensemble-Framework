library(dplyr)
library(MCMCpack)

source('wrapper/pwe_get_trt_effect.R')

## treatment effect (estimand): difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

#' function to estimate the treatment effect (difference in 2-year survival probability between treated and untreated),
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for a CurePWE model under priors other than PSIPP
#'
#' @param post.samples      posterior samples of a CurePWE model under various priors (other than PSIPP), with an attribute
#'                          called 'data' which includes the list of variables specified in the data block of the Stan 
#'                          program.
#' @param trt.name          name of treatment indicator in the data
#'
get.surv.diff.2yr.curepwe <- function(
    post.samples,
    trt.name = "treatment"
) {
  p_cured = as.numeric( 
    suppressWarnings(unlist( post.samples[, 'p_cured'] )) 
  )
  ## compute the difference in 2-year survival probability from a PWE model under priors other than PSIPP
  surv.diff = get.surv.diff.2yr.pwe(
    post.samples = post.samples,
    trt.name = trt.name
  )
  
  return(
    (1 - p_cured) * surv.diff
  )
}


#' function to estimate the treatment effect (difference in 2-year survival probability between treated and untreated),
#' i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0), for a CurePWE model under PSIPP
#'
#' @param post.samples      posterior samples of a PWE model under PSIPP, with an attribute called 'data' which
#'                          includes the list of variables specified in the data block of the Stan program.
#' @param trt.name          name of treatment indicator in the data
#' 
get.surv.diff.2yr.curepwe.psipp <- function(
    post.samples,
    trt.name = "treatment"
) {
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  K         = stan.data$K
  breaks    = stan.data$breaks
  
  # coefficient of treatment indicator for each stratum
  betaMat    = suppressWarnings(
    as.matrix( post.samples[, paste0(trt.name, '_stratum_', 1:K), drop=F] )
  )
  lambdaMat  = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"), drop=F] )
  )
  # cure proportion for each stratum
  p_curedMat = suppressWarnings(
    as.matrix( post.samples[, paste0("p_cured", "_stratum_", 1:K), drop=F] )
  )
  
  ## interval index for which the time t = 2 belongs to
  interval.id = findInterval(2, breaks, left.open = TRUE)
  
  ## compute cumulative baseline hazard at 2 years for each stratum
  cumblhaz = sapply(1:K, function(k){
    lambda = lambdaMat[, paste0("basehaz", "_stratum_", k, "[", 1:J, "]"), drop = F]
    
    ## Compute cumulative baseline hazard at each interval
    H = apply(lambda, 1, function(x){
      as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
    })
    H = matrix(H, nrow = J-1)
    H = cbind(0, t(H))
    return(
      lambda[, interval.id] * (2 - breaks[interval.id]) + H[, interval.id]
    )
  })
  
  ## predicted 2-year survival probability for untreated for each stratum
  S.ctl = lapply(1:K, function(k){
    H_stratum     = cumblhaz[, k]
    S_ctl_stratum = exp( -H_stratum )
    return(S_ctl_stratum)
  }
  )
  S.ctl = do.call(cbind, S.ctl)
  
  ## predicted 2-year survival probability for untreated for each stratum
  S.trt = lapply(1:K, function(k){
    H_stratum     = cumblhaz[, k]
    beta_stratum  = betaMat[, k]
    S_trt_stratum = exp( -H_stratum * exp(beta_stratum) )
    return(S_trt_stratum)
  }
  )
  S.trt = do.call(cbind, S.trt)
  
  surv.diff = (1 - p_curedMat) * (S.trt - S.ctl)
  ## take sample mean across strata
  surv.diff = rowMeans(surv.diff)
  return(surv.diff)
}
