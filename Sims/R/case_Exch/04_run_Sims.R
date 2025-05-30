## Simulation study under exchangeable case
## both true data generating and analysis models are Weibull AFT models
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
wrapper.dir <- "wrapper"
source(file.path(wrapper.dir, "get_strata_data.R"))
source(file.path(wrapper.dir, "sample_ensemble.R"))
source(file.path(wrapper.dir, "aft_get_trt_effect.R"))
source(file.path(wrapper.dir, "generate_data_Exch.R"))

## directory to store files
save.dir   <- '/work/users/x/i/xinxinc/super_prior/Sims/case_Exch'
grid       <- readRDS(file = "Sims/R/case_Exch/grid.rds")
seeds_list <- readRDS(file = "Sims/R/case_Exch/seeds_list.rds")


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
fmla       <- Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
fmla.psipp <- Surv(failtime, failcens) ~ treatment
ps.formula <- ~ sex + cage + node_bin
nStrata    <- 5 ## number of strata for analysis using PSIPP

## setup MCMC parameters
nburnin  <- 2000
nsamples <- 20000
nchains  <- 1

for ( j in seq_len(ndatasets) ) {
  set.seed(seeds.id[j]) # set seed for generating current data
  ## generate current and historical time-to-event data sets
  data.list <- sim.PP(
    prop.cens     = cens.prop.id,
    nevents       = nevents.id,
    n0events      = n0events.id,
    X             = design.mtx.curr,
    theta         = mle.curr,
    trt.prob      = 0.5
  )

  ## analysis using a non-informative reference prior, PP w/ a0 = 0.5, BHM, LEAP, and PSIPP

  ## reference prior:
  d.ref         <- aft.post(
    fmla, data.list
    , dist = "weibull"
    , beta.mean = NULL, beta.sd = NULL
    , scale.mean = NULL, scale.sd = NULL
    , get.loglik = T,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.ref <- aft.logml.post(
    post.samples = d.ref
  )$logml

  ## PP w/ a0 = 0.5:
  d.pp        <- aft.pp(
    fmla, data.list
    , a0 = a0.id
    , dist = "weibull"
    , beta.mean = NULL, beta.sd = NULL
    , scale.mean = NULL, scale.sd = NULL
    , get.loglik = T,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.pp <- aft.logml.pp(
    post.samples = d.pp,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )$logml

  ## BHM:
  d.bhm        <- aft.bhm(
    fmla, data.list
    , dist = "weibull"
    , meta.mean.mean = NULL
    , meta.mean.sd = NULL
    , meta.sd.mean = 0
    , meta.sd.sd = 0.5
    , scale.mean = NULL, scale.sd = NULL
    , get.loglik = T,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.bhm <- aft.logml.map(
    post.samples = d.bhm,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )$logml
  
  ## LEAP:
  d.leap           <- aft.leap(
    fmla, data.list
    , dist = "weibull"
    , K = 2
    , prob.conc = NULL
    , beta.mean = NULL, beta.sd = NULL
    , scale.mean = NULL, scale.sd = NULL
    , gamma.lower = 0
    , gamma.upper = 1
    , get.loglik = T,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.leap   <- aft.logml.leap(
    post.samples = d.leap,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )$logml
  
  ## PSIPP:
  nBorrow    <- nrow(data.list[[2]])
  res.strata <- get.strata.data(
    data            = data.list[[1]], 
    histdata        = data.list[[2]],
    ps_fml_covs     = ps.formula,
    v_arm           = NULL, # treat as a single-arm study since we want to borrow on both arms
    ctl_arm_level   = NULL,
    borrow_ctl_only = FALSE,
    nstrata         = nStrata,
    total_borrow    = nBorrow
  )
  data.list.PSIPP   <- res.strata$data.list
  strata.list       <- res.strata$strata.list

  d.psipp           <- aft.stratified.pp(
    formula = fmla.psipp,
    data.list = data.list.PSIPP
    , strata.list = strata.list
    , a0.strata = res.strata$a0.strata
    , dist = "weibull"
    , beta.mean = NULL, beta.sd = NULL
    , scale.mean = NULL, scale.sd = NULL
    , get.loglik = T,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )
  res.logml.psipp <- aft.logml.stratified.pp(
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

  loo.list  <- lapply(d.list, function(f){
    loglik <- suppressWarnings( f %>%
                                  dplyr::select(starts_with("log_lik")) %>%
                                  as.matrix )
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
          trt.name = "treatment"
        )
      )
    }else{
      return(
        get.surv.diff.2yr(
          post.samples = d,
          trt.name = "treatment"
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
    samples <- sample.ensemble(wts = wts, samples.mtx = trt.effect.mtx)
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

  if(j == 1) {
    simres.trt <- trt.estim.j
    simres.wts <- wts.estim.j
  } else{
    simres.trt <- rbind(simres.trt, trt.estim.j)
    simres.wts <- rbind(simres.wts, wts.estim.j)
  }
  
  print(paste0("#################### Completed iteration ", j, " ######################"))
}

## save the results
lst <- list(
  'id'      = id
  , 'scen'  = grid.id
  , 'simres.trt'  = simres.trt
  , 'simres.wts'  = simres.wts
)
saveRDS(lst, filename)
