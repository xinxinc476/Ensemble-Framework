## Analyze E1690 (current) and E1684 (historical) data sets using a piecewise exponential (PWE) model 
## and a standard cure rate (CurePWE) model with varying number of time intervals/breaks in the PWE 
## model and various priors
library(hdbayes)
library(dplyr)
library(survival)
library(loo)
library(posterior)

## Sampling parameters
iter_warmup   <- 5000
iter_sampling <- 25000
chains        <- 1

## Obtain scenarios for computation
grid <- readRDS("Analysis/R/grid.rds")
## Get scenario for current ID
id <- as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if ( is.na(id) )
  id <- 1

## Directory to save results
save.dir <- 'Analysis/Results'

scen     <- grid[id, ]
prior.id <- scen$priors
model.id <- scen$models
J.id     <- scen$J
## Obtain file name based on id
filename <- file.path(save.dir, 
                      paste0('id_', id, "_", model.id, "_", prior.id, '_nintervals_', J.id, '.rds'))

## load data
hist <- E1684
curr <- E1690

## replace 0 failure times with 0.50 days
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 1 subject w/ failtime = 0
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 10 subjects w/ failtime = 0

## Center and scale age based on current data
age_stats <- with(curr,
                  c('mean' = mean(age), 'sd' = sd(age)))
hist$cage <- ( hist$age - age_stats['mean'] ) / age_stats['sd']
curr$cage <- ( curr$age - age_stats['mean'] ) / age_stats['sd']
data.list <- list(curr, hist)

## formula for current and historical data
fmla       <- survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
fmla.psipp <- survival::Surv(failtime, failcens) ~ treatment

## obtain cut points for intervals
nbreaks <- J.id
probs   <- 1:nbreaks / nbreaks
breaks  <- curr %>%
  filter(failcens == 1) %>%
  reframe(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks  <- c(0, breaks)
breaks[length(breaks)] <- max(10000, 1000 * breaks[length(breaks)])

## source wrapper functions
if( prior.id == "psipp" ){
  wrapper.dir <- 'wrapper'
  source(file.path(wrapper.dir, "get_strata_data.R"))
  ps.formula <- ~ sex + cage + node_bin
  nStrata    <- 4
  nBorrow    <- nrow(hist)
  res.strata <- get.strata.data(
    data = curr, histdata = hist,
    ps_fml_covs     = ps.formula,
    v_arm           = NULL, # treat as a single-arm study since we want to borrow on both arms
    ctl_arm_level   = NULL,
    borrow_ctl_only = FALSE,
    nstrata         = nStrata,
    total_borrow    = nBorrow
  )
  data.list.PSIPP   <- res.strata$data.list
  strata.list       <- res.strata$strata.list
}

if ( model.id == "PWE" ){
  ## fit PWE model
  if ( prior.id == "psipp" ){
    d <- pwe.stratified.pp(
      formula = fmla.psipp,
      data.list = data.list.PSIPP,
      strata.list = strata.list,
      breaks = breaks,
      a0.strata = res.strata$a0.strata,
      beta.mean = 0, beta.sd = 10,
      base.hazard.mean = 0, base.hazard.sd = 10,
      get.loglik = T,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )
    res.logml <- pwe.logml.stratified.pp(
      d,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )
    
  }else {
    if ( prior.id == "ref" ){
      d <- pwe.post(
        formula = fmla,
        data.list = list(curr),
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.post(d)
      
    }else if( prior.id == "pp" ){
      ## PP with a0 = 0.5
      d <- pwe.pp(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        a0 = 0.5,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.pp(
        d,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "bhm" ){
      d <- pwe.bhm(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        meta.mean.mean = 0, meta.mean.sd = 10,
        meta.sd.mean = 0, meta.sd.sd = 0.5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.map(
        d,
        chains = chains, iter_warmup = iter_warmup + 3000, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "cp" ){
      d <- pwe.commensurate(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        beta0.mean = 0, beta0.sd = 10,
        p.spike = 0.1,
        spike.mean = 200, spike.sd = 0.1,
        slab.mean = 0, slab.sd = 5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.commensurate(
        d,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "npp" ){
      a0 <- seq(0, 1, length.out = 101)
      
      library(parallel)
      ncores   <- 5
      ## wrapper to obtain log normalizing constant
      logncfun <- function(a0, ...){
        hdbayes::pwe.npp.lognc(
          formula = fmla, histdata = data.list[[2]], breaks = breaks, a0 = a0,
          beta.mean = 0, beta.sd = 10,
          base.hazard.mean = 0, base.hazard.sd = 10,
          ...
        )
      }
      cl <- makeCluster(ncores)
      clusterSetRNGStream(cl, 123)
      clusterExport(cl, varlist = c('fmla', 'data.list', 'breaks'))
      a0.lognc <- parLapply(
        cl = cl, X = a0, fun = logncfun, 
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
        refresh = 0
      )
      stopCluster(cl)
      a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
      
      d <- pwe.npp(
        formula = fmla,
        data.list = data.list,
        a0.lognc = a0.lognc$a0,
        lognc = a0.lognc$lognc,
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        a0.shape1 = 1, a0.shape2 = 1,
        a0.lower = 0, a0.upper = 1,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.npp(d)
      
    }else if( prior.id == "leap" ){
      d <- pwe.leap(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        K = 2, prob.conc = 1,
        gamma.lower = 0, gamma.upper = 1,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- pwe.logml.leap(
        d,
        chains = chains, iter_warmup = iter_warmup + 3000, iter_sampling = iter_sampling
      )
      
    }else{
      stop("prior.id incorrect")
    }
  }
  
}else{
  ## fit CurePWE model
  if ( prior.id == "psipp" ){
    d <- curepwe.stratified.pp(
      formula = fmla.psipp,
      data.list = data.list.PSIPP,
      strata.list = strata.list,
      breaks = breaks,
      a0.strata = res.strata$a0.strata,
      beta.mean = 0, beta.sd = 10,
      base.hazard.mean = 0, base.hazard.sd = 10,
      logit.pcured.mean = 0, logit.pcured.sd = 3,
      get.loglik = T,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )
    res.logml <- curepwe.logml.stratified.pp(
      d,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )
    
  }else {
    if ( prior.id == "ref" ){
      d <- curepwe.post(
        formula = fmla,
        data.list = list(curr),
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.post(d)
      
    }else if( prior.id == "pp" ){
      ## PP with a0 = 0.5
      d <- curepwe.pp(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        a0 = 0.5,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.pp(
        d,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "bhm" ){
      d <- curepwe.bhm(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        meta.mean.mean = 0, meta.mean.sd = 10,
        meta.sd.mean = 0, meta.sd.sd = 0.5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.map(
        d,
        chains = chains, iter_warmup = iter_warmup + 3000, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "cp" ){
      d <- curepwe.commensurate(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        beta0.mean = 0, beta0.sd = 10,
        p.spike = 0.1,
        spike.mean = 200, spike.sd = 0.1,
        slab.mean = 0, slab.sd = 5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.commensurate(
        d,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "npp" ){
      a0       <- seq(0, 1, length.out = 101)
      
      library(parallel)
      ncores   <- 5
      ## wrapper to obtain log normalizing constant
      logncfun <- function(a0, ...){
        hdbayes::curepwe.npp.lognc(
          formula = fmla, histdata = data.list[[2]], breaks = breaks, a0 = a0,
          beta.mean = 0, beta.sd = 10,
          base.hazard.mean = 0, base.hazard.sd = 10,
          logit.pcured.mean = 0, logit.pcured.sd = 3,
          ...
        )
      }
      cl <- makeCluster(ncores)
      clusterSetRNGStream(cl, 123)
      clusterExport(cl, varlist = c('fmla', 'data.list', 'breaks'))
      a0.lognc <- parLapply(
        cl = cl, X = a0, fun = logncfun, 
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
        refresh = 0
      )
      stopCluster(cl)
      a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
      
      d <- curepwe.npp(
        formula = fmla,
        data.list = data.list,
        a0.lognc = a0.lognc$a0,
        lognc = a0.lognc$lognc,
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        a0.shape1 = 1, a0.shape2 = 1,
        a0.lower = 0, a0.upper = 1,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.npp(d)
      
    }else if( prior.id == "leap" ){
      d <- curepwe.leap(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        K = 2, prob.conc = 1,
        gamma.lower = 0, gamma.upper = 1,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      res.logml <- curepwe.logml.leap(
        d,
        chains = chains, iter_warmup = iter_warmup + 3000, iter_sampling = iter_sampling
      )
      
    }else{
      stop("prior.id incorrect")
    }
  }

}

loglik <- suppressWarnings( 
  d %>% 
    dplyr::select(starts_with("log_lik")) %>% 
    as.matrix 
)
## compute elpd_loo
res.loo <- loo::loo(loglik)

res <- list(
  'scen'      = scen
  , 'draws'     = d
  , 'logml'     = res.logml$logml
  , 'res.loo'   = res.loo
  , 'res.logml' = res.logml
)
if ( prior.id == "psipp" ){
  res <- append(
    res,
    list("res.strata" = res.strata)
  )
}
saveRDS(res, filename)
