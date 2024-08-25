library(hdbayes)
library(cmdstanr)

#' Log marginal likelihood of a GLM under stratified power prior (PP)
#'
#' @param post.samples      output from [glm.stratified.pp()] giving posterior samples of a GLM under the stratified power
#'                          prior (PP), with an attribute called 'data' which includes the list of variables specified in
#'                          the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup`
#'                          in `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method
#'                          in cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
glm.logml.stratified.pp = function(
    post.samples,
    bridge.args   = NULL,
    iter_warmup   = 1000,
    iter_sampling = 1000,
    chains        = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = stan.data$p
  K        = stan.data$K
  X        = stan.data$X
  oldnames = paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  newnames = paste0( colnames(X), '_stratum_', rep(1:K, each = p) )
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, paste0( 'dispersion', '_stratum_', 1:K ))
  }
  colnames(d)[colnames(d) %in% newnames] = oldnames
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion parameters
  stan.data$lognc_disp  = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

  ## compute lognc for stratified PP
  start.idx.hist = stan.data$start_idx_hist
  end.idx.hist   = stan.data$end_idx_hist
  n0.stratum     = end.idx.hist - start.idx.hist + 1 ## sample size for each stratum in historical data
  X0             = stan.data$X0
  y0             = stan.data$y0
  offs0          = stan.data$offs0

  res.lognc.pp = lapply(1:K, function(k){
    ## get Stan data for PP using historical data from each stratum
    y0.stratum    = y0[ start.idx.hist[k]:end.idx.hist[k] ]
    X0.stratum    = X0[ start.idx.hist[k]:end.idx.hist[k], ]
    offs0.stratum = offs0[ start.idx.hist[k]:end.idx.hist[k] ]
    hist.stan.data = list(
      'K'               = 1,
      'N'               = n0.stratum[k],
      'start_idx'       = 1,
      'end_idx'         = n0.stratum[k],
      'p'               = p,
      'y'               = y0.stratum,
      'X'               = X0.stratum,
      'mean_beta'       = stan.data$beta_mean[, k],
      'sd_beta'         = stan.data$beta_sd[, k],
      'a0_vals'         = stan.data$a0s[k],
      'disp_mean'       = stan.data$disp_mean[k],
      'disp_sd'         = stan.data$disp_sd[k],
      'dist'            = stan.data$dist,
      'link'            = stan.data$link,
      'offs'            = offs0.stratum
    )
    ## fit PP using historical data sets
    glm_pp = cmdstanr::cmdstan_model("~/Documents/UNC/Dissertation/super prior/Analysis/Stan/glm_pp.stan")
    fit    = glm_pp$sample(data = hist.stan.data,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                        ...)

    if ( hist.stan.data$dist > 2 ) {
      ## rename parameters
      oldnames = 'dispersion[1]'
      newnames = 'dispersion'
      hist.post.samples = hdbayes:::rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
    }else {
      hist.post.samples = fit$draws(format = 'draws_df')
    }
    attr(x = hist.post.samples, which = 'data') = hist.stan.data

    ## compute log normalizing constant for PP using historical data sets
    res.hist = hdbayes:::glm.pp.lognc(
      post.samples   = hist.post.samples,
      bridge.args    = bridge.args
    )
    return(res.hist)
  })
  lognc_stratified_pp = sum( unlist( lapply(res.lognc.pp, function(l){l[['lognc']]}) ) )

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta       = pars[paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta       = matrix(beta, nrow = p, ncol = K)
    prior_lp   = sum( sapply(1:K, function(k){
      dnorm(beta[, k], mean = data$beta_mean[, k], sd = data$beta_sd[, k], log = T)
      })
    )
    dist       = data$dist
    link       = data$link
    dispersion = rep(1.0, K)
    if ( dist > 2 ){
      dispersion = pars[paste0( 'dispersion[', 1:K, ']' )]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
    }

    data_lp    = sum( sapply(1:K, function(k){
      y         = data$y[  data$start_idx_curr[k]:data$end_idx_curr[k] ]
      X         = data$X[ data$start_idx_curr[k]:data$end_idx_curr[k], ]
      offs      = data$offs[ data$start_idx_curr[k]:data$end_idx_curr[k] ]
      y0        = data$y0[  data$start_idx_hist[k]:data$end_idx_hist[k] ]
      X0        = data$X0[ data$start_idx_hist[k]:data$end_idx_hist[k], ]
      offs0     = data$offs0[ data$start_idx_hist[k]:data$end_idx_hist[k] ]

      hdbayes:::glm_lp(y, beta[, k], X, dist, link, offs, dispersion[k]) +
        data$a0s[k] * hdbayes:::glm_lp(y0, beta[, k], X0, dist, link, offs0, dispersion[k])
      })
    )
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p*K)
  ub = rep(Inf, p*K)
  if( stan.data$dist > 2 ) {
    lb = c(lb, rep(0, K))
    ub = c(ub, rep(Inf, K))
  }
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

  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## and outputs from computing the log normalizing constant for stratified PP
  res = list(
    'model'        = "stratified PP",
    'logml'        = bs$logml - lognc_stratified_pp,
    'bs'           = bs,
    'res.lognc.pp' = res.lognc.pp
  )
  return(res)
}
