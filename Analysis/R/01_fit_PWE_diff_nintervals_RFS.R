## Fit a proportional hazards (PH) model with a piece-wise constant baseline hazards using the equivalent Poisson likelihood
## using different number of time intervals/breaks and various priors
remove(list = ls())

library(tidyverse)
library(survival)
library(hdbayes)
library(loo)
library(posterior)
library(parallel)

## Sampling parameters
nburnin  <- 5000
nsamples <- 25000
nchains  <- 1

## Obtain scenarios for computation
priors    <- c('ref', 'pp_0.2', 'pp_0.5', 'pp_0.8',  'psipp', 'bhm', 'commensurate', 'leap', 'npp')
J         <- 1:10
scenarios <- expand.grid(
  J = J, priors = priors
  , stringsAsFactors = FALSE
)

## Get scenario for current ID
id <- as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if ( is.na(id) )
  id <- 1

## Directory to save results
save.dir <- '/work/users/x/i/xinxinc/super_prior/PWE_choose_nintervals_RFS'

scen     <- scenarios[id, ]
prior.id <- scen$priors
J.id     <- scen$J

## Obtain file name based on id
filename <- file.path(save.dir, paste0('PWE_id_', id, "_", prior.id, '_nintervals_', J.id, '.rds'))


## source wrappers
wrapper.dir <- '~/projects/super_prior/Analysis/R/wrapper_glm'
source(file.path(wrapper.dir, 'wrapper_surv.R'))
source(file.path(wrapper.dir, 'get_loglik.R'))
if (prior.id == "psipp"){
  source(file.path(wrapper.dir, 'get_strata_data.R'))
  source(file.path(wrapper.dir, 'glm_stratified_pp.R'))
  source(file.path(wrapper.dir, 'glm_logml_stratified_pp.R'))
}

## load data
hist <- E1684
curr <- E1690

## replace 0 failure times with 0.50 days
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 1 subject w/ failtime = 0
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 10 subjects w/ failtime = 0

## replace 0 overall survival times with 0.50
#hist <- hist %>% mutate(survtime = if_else(survtime == 0, 0.50, survtime)) # 1 subject w/ survtime = 0
#no subjects in current data w/ survtime = 0

## Center and scale age based on current data
age_stats <- with(curr,
                  c('mean' = mean(age), 'sd' = sd(age)))
hist$cage <- ( hist$age - age_stats['mean'] ) / age_stats['sd']
curr$cage <- ( curr$age - age_stats['mean'] ) / age_stats['sd']

#hist$cage <- as.numeric( scale(hist$age) )
#curr$cage <- as.numeric( scale(curr$age) )

## obtain cut points for intervals
nbreaks <- J.id
probs   <- 1:nbreaks / nbreaks
breaks  <- curr %>%
  filter(failcens == 1) %>%
  reframe(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks  <- breaks[-J.id]

## formula for current and historical data (including intercept)
fmla       <- failcens ~ interval + treatment + sex + cage + node_bin
fmla.psipp <- failcens ~ interval + treatment

if (J.id == 1 ){
  fmla       <- failcens ~ treatment + sex + cage + node_bin
  fmla.psipp <- failcens ~ treatment
}
family <- poisson("log")

if ( prior.id == "psipp" ){
  ps.formula <- ~ sex + cage + node_bin
  nStrata    <- 4
  nBorrow    <- nrow(hist)
  res.strata <- get.strata.data(
    data.list       = list(curr, hist),
    ps_fml_covs     = ps.formula,
    v_arm           = "treatment",
    ctl_arm_level   = 0,
    borrow_ctl_only = FALSE,
    nstrata         = nStrata,
    total_borrow    = nBorrow
  )
  data.list   <- res.strata$data.list
  strata.list <- res.strata$strata.list
  data.list   <- lapply(1:length(data.list), function(i){
    df = data.list[[i]]
    df$ps_strata = strata.list[[i]]
    return(df)
  })

  ## transform data into counting process format
  covariate.name <- c("treatment", "sex", "cage", "node_bin", "ps_strata")
  curr.poisson   <- get.counting.process.data(
    data = data.list[[1]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
    covariate.name = covariate.name
  )
  hist.poisson   <- get.counting.process.data(
    data = data.list[[2]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
    covariate.name = covariate.name
  )
  data.list      <- list(currdata = curr.poisson, histdata = hist.poisson)
  offset.list    <- list(offs = curr.poisson$log_exposue, offs0 = hist.poisson$log_exposue)
  strata.list    <- lapply(data.list, function(l){
    l[, 'ps_strata']
  })

  d <- glm.stratified.pp(
    formula = fmla.psipp, family = family, data.list = data.list,
    strata.list = strata.list,
    a0.strata = res.strata$a0.strata,
    offset.list = offset.list,
    is.sorted = FALSE,
    iter_warmup = nburnin,
    iter_sampling = nsamples,
    chains = nchains
  )

  loglik    <- get.loglik.glm.stratified.pp(post.samples = d)
  res.logml <- glm.logml.stratified.pp(
    post.samples = d,
    iter_warmup = nburnin + 3000,
    iter_sampling = nsamples + 3000,
    chains = nchains
  )

}else {
  ## transform data into counting process format
  covariate.name <- c("treatment", "sex", "cage", "node_bin")
  curr.poisson   <- get.counting.process.data(
    data = curr, breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
    covariate.name = covariate.name
  )
  hist.poisson   <- get.counting.process.data(
    data = hist, breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
    covariate.name = covariate.name
  )
  data.list      <- list(currdata = curr.poisson, histdata = hist.poisson)
  offset.list    <- list(offs = curr.poisson$log_exposue, offs0 = hist.poisson$log_exposue)

  if ( prior.id == "ref" ){
    d <- glm.post(
      formula = fmla, family = family, data.list = data.list,
      offset.list = offset.list,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.post(
      post.samples = d
    )

  }else if( prior.id == "pp_0.2" ){
    d <- glm.pp(
      formula = fmla, family = family, data.list = data.list,
      a0.vals = 0.2, offset.list = offset.list,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.pp(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains
    )

  }else if( prior.id == "pp_0.5" ){
    d <- glm.pp(
      formula = fmla, family = family, data.list = data.list,
      a0.vals = 0.5, offset.list = offset.list,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.pp(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains
    )

  }else if( prior.id == "pp_0.8" ){
    d <- glm.pp(
      formula = fmla, family = family, data.list = data.list,
      a0.vals = 0.8, offset.list = offset.list,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.pp(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains
    )

  }else if( prior.id == "bhm" ){
    d <- glm.bhm(
      formula = fmla, family = family, data.list = data.list,
      offset.list = offset.list,
      meta.sd.sd = 0.5,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.map(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains
    )

  }else if( prior.id == "commensurate" ){
    d <- glm.commensurate(
      formula = fmla, family = family, data.list = data.list,
      offset.list = offset.list,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.commensurate(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains
    )

  }else if( prior.id == "npp" ){
    a0       <- seq(0, 1, length.out = 101)
    histdata <- data.list[[2]]

    ## wrapper to obtain log normalizing constant in parallel package
    logncfun <- function(a0, ...){
      hdbayes::glm.npp.lognc(
        formula = fmla, family = family, histdata = histdata, offset0 = histdata$log_exposue, a0 = a0, ...
      )
    }

    #cl <- makeCluster(10)
    #clusterSetRNGStream(cl, 123)
    #clusterExport(cl, varlist = c('fmla', 'family', 'histdata'))
    ## call created function
    #a0.lognc <- parLapply(
    #cl = cl, X = a0, fun = logncfun, iter_warmup = nburnin,
    #iter_sampling = nsamples, chains = nchains, refresh = 0
    #)
    #stopCluster(cl)
    a0.lognc <- lapply(a0, function(a){
      logncfun(a0 = a, iter_warmup = nburnin,
               iter_sampling = nsamples, chains = nchains, refresh = 0)
    })
    a0.lognc <- data.frame( do.call(rbind, a0.lognc) )

    d <- glm.npp(
      formula = fmla, family = family, data.list = data.list,
      offset.list = offset.list,
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml <- glm.logml.npp(
      post.samples = d
    )

  }else if( prior.id == "leap" ){
    K           <- 2
    offset.list <- lapply(offset.list, function(l){
      matrix(rep(l, K), ncol = K)
    })
    d           <- glm.leap(
      formula = fmla, family = family, data.list = data.list,
      offset.list = offset.list,
      K = K,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains
    )
    res.logml   <- glm.logml.leap(
      post.samples = d,
      iter_warmup = nburnin + 3000,
      iter_sampling = nsamples + 3000,
      chains = nchains,
      bridge.args = list(method = "warp3")
    )

  }else{
    stop("prior.id incorrect")
  }

  loglik <- get.loglik.glm(post.samples = d)
}

## compute DIC
dic <- compute.DIC(post.samples = d, is.stratified.pp = (prior.id == "psipp"))

## compute elpd_loo
res.loo <- loo::loo(loglik)

res <- list(
  'draws'       = d
  , 'dic'       = dic
  , 'looic'     = res.loo$estimates['looic', 1]
  , 'logml'     = res.logml$logml
  , 'scen'      = scenarios[id, ]
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
