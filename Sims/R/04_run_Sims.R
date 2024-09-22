library(tidyverse)
library(survival)
library(hdbayes)
library(posterior)
library(cmdstanr)
library(psrwe)
library(bridgesampling)
library(collapse)
library(loo)
library(scoringutils)

## source wrappers
wrapper.sim.dir <- 'R/wrapper_sim'
source(file.path(wrapper.sim.dir, "generate_data.R"))
source(file.path(wrapper.sim.dir, "get_trt_effect.R"))

wrapper.glm.dir <- 'R/wrapper_glm'
source(file.path(wrapper.glm.dir, 'glm_stratified_pp.R'))
source(file.path(wrapper.glm.dir, 'glm_logml_stratified_pp.R'))

## directory to store files
save.dir   <- 'Sims/Results'
grid       <- readRDS(file = "Sims/R/grid.rds")
seeds_list <- readRDS(file = "Sims/R/seeds_list.rds")


## get task ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id ) )
  id <- 1

## get this job's sim parameters
grid.id      <- grid[id, ]
seeds.id     <- seeds_list[grid.id$start:grid.id$end]
rm(grid, seeds_list)
nevents.id   <- as.integer(grid.id$nevents)
n0events.id  <- as.integer(grid.id$n0events)
cens.prop.id <- as.numeric(grid.id$cens.prop)
a0.id        <- 0.5
case.id      <- as.character(grid.id$case)
ndatasets    <- grid.id$end - grid.id$start + 1


## obtain file name based on id
filename     <- file.path(save.dir, paste0('id_', id, '_', 'nevents_', nevents.id,
                                           '_censorProp_', cens.prop.id,
                                           '_case_', case.id,
                                           '_rep_', grid.id$start, '-', grid.id$end, '.rds'))

## setup simulation parameters
nbreaks    <- 5
## formula for current and historical data (including intercept)
fmla       <- failcens ~ interval + treatment + sex + cage + node_bin
fmla.psipp <- failcens ~ interval + treatment
family     <- poisson("log")
ps.formula <- ~ sex + cage + node_bin
nStrata    <- 3 ## number of strata for analysis using PSIPP

## setup MCMC parameters
nburnin  <- 5000
nsamples <- 25000
nchains  <- 1


for ( j in seq_len(ndatasets) ) {
  set.seed(seeds.id[j]) # set seed for generating current data

  ## generate current and historical time-to-event data sets
  if( case.id == "UNEXCH" ){
    data.list <- sim.UNEXCH(
      prop.cens     = cens.prop.id,
      nevents       = nevents.id,
      n0events      = n0events.id,
      X             = design.mtx.curr,
      X0            = design.mtx.hist,
      theta         = mle.curr,
      theta0        = mle.hist,
      trt.prob      = 0.5
    )
  }else if( case.id == "PP" ){
    data.list <- sim.PP(
      prop.cens     = cens.prop.id,
      nevents       = nevents.id,
      n0events      = n0events.id,
      X             = design.mtx.curr,
      theta         = mle.curr,
      trt.prob      = 0.5
    )
  }else if( case.id == "BHM" ){
    data.list <- sim.BHM(
      prop.cens     = cens.prop.id,
      nevents       = nevents.id,
      n0events      = n0events.id,
      X             = design.mtx.curr,
      theta         = mle.curr,
      theta.sd      = 0.2,
      trt.prob      = 0.5
    )
  }else if( case.id == "LEAP" ){
    data.list <- sim.LEAP(
      prop.cens     = cens.prop.id,
      nevents       = nevents.id,
      n0events      = n0events.id,
      X             = design.mtx.curr,
      X0            = design.mtx.hist,
      theta         = mle.curr,
      theta0        = mle.hist,
      exch          = 0.5,
      trt.prob      = 0.5
    )
  }else if( case.id == "PSIPP" ){
    data.list <- sim.PSIPP(
      prop.cens     = cens.prop.id,
      nevents       = nevents.id,
      n0events      = n0events.id,
      X             = design.mtx.curr,
      X0            = design.mtx.hist,
      ps.formula    = ~ sex + cage + node_bin,
      nStrata       = 3,
      theta         = mle.curr,
      drift.beta    = c(-0.2, 0, 0.2),
      trt.prob      = 0.5
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

  ## analysis using a non-informative reference prior, PP w/ a0 = 0.5, BHM, LEAP, and PSIPP

  ## reference prior:
  d.ref         <- glm.post(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = input.list.nonpsipp$offset.list,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.ref <- glm.logml.post(
    post.samples = d.ref
  )$logml

  ## PP w/ a0 = 0.5:
  d.pp        <- glm.pp(
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
  )$logml

  ## BHM:
  d.bhm        <- glm.bhm(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = input.list.nonpsipp$offset.list,
    meta.sd.sd = 0.5,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.bhm <- glm.logml.map(
    post.samples = d.bhm,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )$logml

  ## LEAP:
  K                <- 2
  offset.list.leap <- lapply(input.list.nonpsipp$offset.list, function(l){
    matrix(rep(l, K), ncol = K)
  })
  d.leap           <- glm.leap(
    formula = fmla, family = family, data.list = input.list.nonpsipp$data.list,
    offset.list = offset.list.leap,
    K = K,
    gamma.upper = 1,
    gamma.lower = 0,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.leap   <- glm.logml.leap(
    post.samples = d.leap,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )$logml

  ## PSIPP:
  d.psipp        <- glm.stratified.pp(
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
  )$logml

  d.list <- list(
    "ref" = d.ref, "pp" = d.pp, "bhm" = d.bhm,
    "leap" = d.leap, "psipp" = d.psipp
  )
  rm(d.ref, d.bhm, d.leap, d.pp, d.psipp)

  ## estimated weights from each model-averaging method
  wts.estim.j <- data.frame(
    method  = names(d.list)
  )
  wts.estim.j$logml     <- c(
    res.logml.ref, res.logml.pp, res.logml.bhm,
    res.logml.leap, res.logml.psipp
  )
  wts.estim.j$post.prob <- bridgesampling::post_prob(wts.estim.j$logml)

  loo.list  <- lapply(1:length(d.list), function(i){
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
  lpd_point <- sapply(loo.list, function(l){
    l$pointwise[,"elpd_loo"]
  })
  wts.estim.j$pbma_wts     <- pseudobma_weights(lpd_point, BB = FALSE)
  wts.estim.j$pbma_BB_wts  <- pseudobma_weights(lpd_point, BB = TRUE) # default is BB=TRUE
  wts.estim.j$stacking_wts <- stacking_weights(lpd_point)
  rm(loo.list, lpd_point)


  ## obtain posterior samples of treatment effects from each individual prior
  trt.effect.list <- lapply(1:length(d.list), function(i){
    d     <- d.list[[i]]
    prior <- names(d.list)[i]
    if( prior == "psipp" ){
      return(
        get.surv.diff.2yr.psipp(
          post.samples = d,
          breaks = input.list.psipp$breaks
        )
      )
    }else{
      return(
        get.surv.diff.2yr.glm(
          post.samples = d,
          data = data.list$curr,
          breaks = input.list.nonpsipp$breaks
        )
      )
    }
  })
  trt.estim.j <- sapply(1:length(trt.effect.list), function(i){
    s <- trt.effect.list[[i]]
    c(mean = mean(s), sd = sd(s), posterior::quantile2(s, probs = c(0.025, 0.975)))
  })
  trt.estim.j        <- t(trt.estim.j) %>%
    as.data.frame()
  trt.estim.j$method <- names(d.list)
  trt.estim.j        <- trt.estim.j[, c("method", colnames(trt.estim.j)[which(colnames(trt.estim.j) != "method")])]

  ## obtain posterior samples of treatment effects from each model-averaged prior
  trt.effect.mtx     <- do.call(cbind, trt.effect.list)
  trt.effect.list    <- lapply(c("post.prob", "pbma_wts", "pbma_BB_wts", "stacking_wts"), function(wts_name){
    wts     <- as.numeric( wts.estim.j[, wts_name] )
    samples <- sample.model.avg.prior(wts = wts, samples.mtx = trt.effect.mtx)
    return(samples)
  })
  trt.estim.j2 <- sapply(1:length(trt.effect.list), function(i){
    s <- trt.effect.list[[i]]
    c(mean = mean(s), sd = sd(s), posterior::quantile2(s, probs = c(0.025, 0.975)))
  })
  trt.estim.j2 <- t(trt.estim.j2) %>%
    as.data.frame()
  trt.estim.j2$method <- c("bma", "pbma", "pbma_BB", "stacking")
  trt.estim.j2        <- trt.estim.j2[, c("method", colnames(trt.estim.j2)[which(colnames(trt.estim.j2) != "method")])]
  trt.estim.j         <- rbind(trt.estim.j, trt.estim.j2) %>%
    as.data.frame()
  rm(trt.estim.j2, trt.effect.list, trt.effect.mtx)

  ## store results
  if(j == 1) {
    simres.trt <- trt.estim.j
    simres.wts <- wts.estim.j
  } else{
    simres.trt <- rbind(simres.trt, trt.estim.j)
    simres.wts <- rbind(simres.wts, wts.estim.j)
  }

  lst <- list(
    'id'      = id
    , 'a0'    = a0.id
    , 'scen'  = grid.id
    , 'start' = grid.id$start
    , 'end' = grid.id$start + j - 1
    , 'simres.trt'  = simres.trt
    , 'simres.wts'  = simres.wts
  )
  temp.dir <- "Sims/Results/in_progress"
  temp.filename <- file.path(temp.dir, paste0('id_', id, '_', 'nevents_', nevents.id,
                                              '_censorProp_', cens.prop.id,
                                              '_case_', case.id,
                                              '_rep_', grid.id$start, '-', grid.id$end, '.rds'))
  saveRDS(lst, temp.filename)
  print(paste0("#################### Completed iteration ", j, " ######################"))

}
## SAVE THE RESULTS
lst <- list(
  'id'      = id
  , 'a0'    = a0.id
  , 'scen'  = grid.id
  , 'simres.trt'  = simres.trt
  , 'simres.wts'  = simres.wts
)
saveRDS(lst, filename)
