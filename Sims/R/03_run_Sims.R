remove(list = ls())

library(tidyverse)
library(survival)
library(hdbayes)
library(posterior)
library(cmdstanr)
library(psrwe)
library(bridgesampling)
library(loo)
library(scoringutils)

## source wrappers
source("Sims/R/wrapper_sim.R")

wrapper.dir <- 'Analysis/R/wrapper_glm'
source(file.path(wrapper.dir, 'wrapper_surv.R'))
source(file.path(wrapper.dir, 'get_loglik.R'))
source(file.path(wrapper.dir, 'get_strata_data.R'))
source(file.path(wrapper.dir, 'glm_stratified_pp.R'))
source(file.path(wrapper.dir, 'glm_logml_stratified_pp.R'))

## directory to store files
save.dir   <- 'Sims/results'
grid       <- readRDS(file = "Sims/R/grid.rds")
seeds_list <- readRDS(file = 'Sims/R/seeds_list.rds')


## get task ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id ) )
  id <- 1

## get this job's sim parameters
grid.id   <- grid[id, ]
rm(grid)
seeds.id  <- seeds_list[grid.id$start:grid.id$end]
n.id      <- as.integer(grid.id$n)
n0.id     <- 250
a0.id     <- 0.5
case.id   <- as.character(grid.id$case)
ndatasets <- grid.id$end - grid.id$start + 1


## obtain file name based on id
filename  <- file.path(save.dir, paste0('id_', id, '_', 'n_', n.id, '_case_', case.id,
                                        '_rep_', grid.id$start, '-', grid.id$end, '.rds'))

## setup simulation parameters
nbreaks    <- 5
## formula for current and historical data (including intercept)
fmla       <- failcens ~ interval + treatment + sex + cage + node_bin
fmla.psipp <- failcens ~ interval + treatment
family     <- poisson("log")
ps.formula <- ~ sex + cage + node_bin
nStrata    <- 4 ## number of strata

## specify true value of treatment effect (log hazard ratio at 2 years for treated v.s. untreated)
param_true <- as.numeric( -mle.curr[3]/exp(mle.curr[1]) )

## setup MCMC parameters
nburnin  <- 2000
nsamples <- 20000
nchains  <- 1


for ( j in seq_len(ndatasets) ) {
  set.seed(seeds.id[j]) # set seed for generating current data
  ## generate data set pair
  if( case.id == "non-exchangeable" ){
    data.list <- sim.UNEXCH(
      n = n.id, n0 = n0.id
    )
  }else if( case.id == "PP" ){
    data.list <- sim.PP(
      n = n.id, n0 = n0.id
    )
  }else if( case.id == "BHM" ){
    data.list <- sim.BHM(
      n = n.id, n0 = n0.id
    )
  }else if( case.id == "LEAP" ){
    data.list <- sim.LEAP(
      n = n.id, n0 = n0.id
    )
  }else{
    stop("not available")
  }

  ## preprocess data.list to be used for fitting different priors
  input.list.psipp <- preprocess.data.list(
    data.list = data.list,
    is.psipp = TRUE,
    nbreaks = nbreaks,
    ps.formula = ps.formula,
    nStrata = nStrata
  )

  input.list.nonpsipp <- preprocess.data.list(
    data.list = data.list,
    is.psipp = FALSE,
    nbreaks = nbreaks,
    ps.formula = ps.formula,
    nStrata = nStrata
  )

  ## analysis using the reference prior, PP w/ a0 = 0.5, BHM, LEAP, and PSIPP

  ## reference prior:
  d.ref <- glm.post(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = input.list.nonpsipp$offset.list,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.ref <- glm.logml.post(
    post.samples = d.ref
  )

  ## PP w/ a0 = 0.5:
  d.pp <- glm.pp(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    a0.vals = a0.id, offset.list = input.list.nonpsipp$offset.list,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.pp <- glm.logml.pp(
    post.samples = d.pp,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )

  ## BHM:
  d.bhm <- glm.bhm(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = input.list.nonpsipp$offset.list,
    meta.sd.sd = 0.5,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.bhm <- glm.logml.map(
    post.samples = d.bhm,
    iter_warmup = nburnin + 3000,
    iter_sampling = nsamples + 3000,
    chains = nchains
  )

  ## LEAP:
  K                <- 2
  offset.list.leap <- lapply(input.list.nonpsipp$offset.list, function(l){
    matrix(rep(l, K), ncol = K)
  })
  d.leap           <- glm.leap(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = offset.list.leap,
    K = K,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.leap   <- glm.logml.leap(
    post.samples = d.leap,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )

  ## PSIPP:
  d.psipp <- glm.stratified.pp(
    formula = fmla.psipp, family = family, data.list = input.list.psipp$data.list,
    strata.list = input.list.psipp$strata.list,
    a0.strata = input.list.psipp$a0.strata,
    offset.list = input.list.psipp$offset.list,
    is.sorted = FALSE,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.psipp <- glm.logml.stratified.pp(
    post.samples = d.psipp,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )

  ## estimated treatment effect from each individual prior
  d.list <- list(
    "ref" = d.ref, "pp" = d.pp, "bhm" = d.bhm,
    "leap" = d.leap, "psipp" = d.psipp
  )
  trt.effect.list <- lapply(1:length(d.list), function(i){
    d     <- d.list[[i]]
    prior <- names(d.list)[i]
    if( prior == "psipp" ){
      return(
        get.log.hazard.ratio.2yr.psipp(
          post.samples = d,
          res.strata = input.list.psipp$res.strata,
          breaks = input.list.psipp$breaks
        )
      )
    }else{
      return(
        get.log.hazard.ratio.2yr.glm(
          post.samples = d,
          data = data.list$curr,
          breaks = input.list.nonpsipp$breaks
        )
      )
    }
  })

  loo.list       <- lapply(1:length(d.list), function(i){
    d     <- d.list[[i]]
    prior <- names(d.list)[i]
    if( prior == "psipp" ){
      loglik <- get.loglik.glm.stratified.pp(post.samples = d)
    }else{
      loglik <- get.loglik.glm(post.samples = d)
    }
    return(
      loo::loo(loglik)
    )
  })
  logml.list     <- list(length = length(file.list))

  trt.estim <- rbind(
    d.ref %>%
      summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975))) %>%
      filter(variable == "trt"),
    d.pp %>%
      summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975))) %>%
      filter(variable == "trt"),
    res$fit$summary("trteff", mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))

  )
  trt.estim$variable     <- c("Ref", "PP", "PSIPP")
  colnames(trt.estim)[1] <- "method"

  ## estimated weights from each method
  wts.estim <- data.frame(
    method  = c("Ref", "PP", "PSIPP")
  )
  wts.estim$logml        <- c(res.logml.ref$logml, res.logml.pp$logml, res.logml.psipp$logml)
  wts.estim$post.prob    <- bridgesampling::post_prob(res.logml.ref$logml, res.logml.pp$logml, res.logml.psipp$logml)

  loo_psipp  <- (res$fit)$loo(k_threshold=0.7)
  ## compute pointwise log-likelihood matrices for PP and Ref
  loo_pp     <- get.Loglik.ref(d.pp, formula, family, data.list)
  loo_ref    <- get.Loglik.ref(d.ref, formula, family, data.list)
  loo_list   <- list(loo_ref, loo_pp, loo_psipp)
  lpd_point  <- sapply(loo_list, function(l){
    l$pointwise[,"elpd_loo"]
  })

  wts.estim$pbma_wts     <- pseudobma_weights(lpd_point, BB=FALSE)
  wts.estim$pbma_BB_wts  <- pseudobma_weights(lpd_point) # default is BB=TRUE
  wts.estim$stacking_wts <- stacking_weights(lpd_point)

  ## obtain posterior samples of treatment effect from each method
  trt.samples.ref = suppressWarnings(
    as.numeric(unlist(d.ref[, "trt"]))
  )
  trt.samples.pp = suppressWarnings(
    as.numeric(unlist(d.pp[, "trt"]))
  )
  trt.samples.psipp <- as.numeric(res$fit$draws(variables = "trteff"))
  trt.samples.mtx   <- cbind(trt.samples.ref, trt.samples.pp, trt.samples.psipp)

  ## draw i.i.d. samples (c0) from categorical distribution with probability being posterior model probabilities
  trt.samples.mtx2 <- sapply(c("post.prob", "pbma_wts", "pbma_BB_wts", "stacking_wts"), function(wts_name){
    wts = as.numeric( wts.estim[, wts_name] )
    samples = replicate(
      iter_sampling,
      sample.BMA(samples.mtx = trt.samples.mtx, post.model.probs = wts)
    )
    return(samples)
  })
  colnames(trt.samples.mtx2)[1] <- "bma"
  trt.estim2 <- apply(trt.samples.mtx2, 2, function(s){
    c(mean = mean(s), sd = sd(s), quantile2(s, probs = c(0.025, 0.975)))
  })
  trt.estim2 <- data.frame(
    method = colnames(trt.samples.mtx2),
    t(trt.estim2),
    row.names = NULL
  )
  rm(trt.samples.mtx, trt.samples.mtx2)
  trt.estim  <- rbind(trt.estim, trt.estim2)

  trt.estim$param_true <- param_true
  ## compute interval score
  ## interval score = 0.05/2 * ( (q97.5-q2.5) + 2/0.05 * max(0, q2.5 - param_true) + 2/0.05 * max(0, param_true - q97.5) )
  trt.estim = trt.estim %>%
    mutate(
      interval_score = scoringutils::interval_score(
        true_values = param_true, lower = `q2.5`, upper = `q97.5`,
        interval_range = 95, weigh = T
      ),
      diff = mean - param_true,
      diff_sq = diff^2,
      ci.ind  = ifelse(param_true >= `q2.5` & param_true <= `q97.5`, 1, 0),
      log_var = ifelse(sd > 0, log(sd^2), NA)
    )

  ## add simulation scenario
  trt.estim$n    <- n.id
  trt.estim$a0   <- a0.id
  trt.estim$case <- case.id
  wts.estim$n    <- n.id
  wts.estim$a0   <- a0.id
  wts.estim$case <- case.id

  ## store results
  if(j == 1) {
    simres.trt <- trt.estim
    simres.wts <- wts.estim
  } else{
    simres.trt <- rbind(simres.trt, trt.estim)
    simres.wts <- rbind(simres.wts, wts.estim)
  }

  if(j %% 10 == 0){
    lst <- list(
      'id'      = id
      , 'scen'  = grid.id[1:3]
      , 'start' = grid.id$start
      , 'end' = grid.id$start + j - 1
      , 'simres.trt'  = simres.trt
      , 'simres.wts'  = simres.wts
    )
    temp.dir <- "/work/users/x/i/xinxinc/super_prior/Sims_small/in_progress"
    temp.filename <- file.path(temp.dir, paste0('id_', id, '_', 'n_', n.id, '_a0_', a0.id, '_case_', case.id,
                                                '_rep_', grid.id$start, '-', grid.id$end, '.rds'))
    saveRDS(lst, temp.filename)
  }
  print(paste0("#################### Completed iteration ", j, " ######################"))

}
## SAVE THE RESULTS
lst <- list(
  'id'      = id
  , 'scen'  = grid.id[1:3]
  , 'simres.trt'  = simres.trt
  , 'simres.wts'  = simres.wts
)
saveRDS(lst, filename)
