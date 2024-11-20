options("mc.cores" = 10)

library(hdbayes)
library(survival)
library(tidyverse)
library(loo)
library(bridgesampling)

wrapper.dir <- "R/wrapper_aft"
source(file.path(wrapper.dir, "wrapper_psipp.R"))
source(file.path(wrapper.dir, "aft_get_trt_effect.R"))
source(file.path(wrapper.dir, "generate_data_PP.R"))

prop.cens <- 1-mean(E1690$failcens)
nevents   <- 200
n0events  <- 200

# simulate survival data under PP assumption
set.seed(123)
data.list <- sim.PP(
  prop.cens     = prop.cens,
  nevents       = nevents,
  n0events      = n0events,
  X             = design.mtx.curr,
  theta         = mle.curr,
  trt.prob      = 0.5
)
data.list <- lapply(1:length(data.list), function(i){
  df    <- data.list[[i]]
  df$ID <- 1:nrow(df)
  return(df)
})

# compute true values of beta and scale (coded in Stan file)
param.true <- mle.curr

fmla     <- Surv(failtime, failcens) ~ treatment + sex + cage + node_bin

warmup  <- 2000
iter    <- 20000 + warmup
chains  <- 1

# reference prior
fit.ref <- aft.post(
  fmla, data.list
  , dist = "weibull"
  , beta.mean = NULL, beta.sd = NULL
  , scale.mean = NULL, scale.sd = NULL
  , get.loglik = T,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
pars.name <- c("scale", "(Intercept)", "treatment", "sex", "cage", "node_bin")
res.ref <- fit.ref  %>%
  dplyr::select(all_of(pars.name)) %>% 
  summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
res.ref$param_true <- param.true
res.ref

# compute log marginal likelihood
logml.ref <- aft.logml.post(
  post.samples = fit.ref
)
summary(logml.ref$bs)
logml.ref.val <- logml.ref$logml

# compare with the corresponding frequentist estimates
fit <- survreg(formula = fmla, data = data.list[[1]], dist = "weibull")
c(fit$scale, fit$coefficients)


# power prior w/ a0 = 0.5
fit.pp <- aft.pp(
  fmla, data.list
  , a0 = 0.5
  , dist = "weibull"
  , beta.mean = NULL, beta.sd = NULL
  , scale.mean = NULL, scale.sd = NULL
  , get.loglik = T,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
res.pp <- fit.pp  %>%
  dplyr::select(all_of(pars.name)) %>% 
  summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
res.pp$param_true <- param.true
res.pp

# compute log marginal likelihood
logml.pp <- aft.logml.pp(
  post.samples = fit.pp
  , bridge.args = NULL
  , iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
summary(logml.pp$bs)
summary(logml.pp$bs.hist)
logml.pp.val <- logml.pp$logml


# BHM
fit.bhm <- aft.bhm(
  fmla, data.list
  , dist = "weibull"
  , meta.mean.mean = NULL
  , meta.mean.sd = NULL
  , meta.sd.mean = 0
  , meta.sd.sd = 0.5
  , scale.mean = NULL, scale.sd = NULL
  , get.loglik = T,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
res.bhm <- fit.bhm %>%
  dplyr::select(all_of(pars.name)) %>% 
  summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
res.bhm$param_true <- param.true
res.bhm

# compute log marginal likelihood
logml.bhm <- aft.logml.map(
  post.samples = fit.bhm
  , bridge.args = NULL
  , iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
summary(logml.bhm$bs)
summary(logml.bhm$bs.hist)
logml.bhm.val <- logml.bhm$logml


# LEAP
fit.leap <- aft.leap(
  fmla, data.list
  , dist = "weibull"
  , K = 2
  , prob.conc = NULL
  , beta.mean = NULL, beta.sd = NULL
  , scale.mean = NULL, scale.sd = NULL
  , gamma.lower = 0
  , gamma.upper = 1
  , get.loglik = T,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
) 
res.leap <- fit.leap %>%
  dplyr::select(all_of(c(pars.name, "probs[1]"))) %>% 
  summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
res.leap$param_true <- c(param.true, 1)
res.leap

# compute log marginal likelihood
logml.leap <- aft.logml.leap(
  post.samples = fit.leap
  , bridge.args = NULL
  , iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
summary(logml.leap$bs)
summary(logml.leap$bs.hist)
logml.leap.val <- logml.leap$logml


# fit PSIPP
nBorrow    <- nrow(data.list[[2]])
nStrata    <- 5
ps.formula <- ~ sex + cage + node_bin
fmla.psipp <- Surv(failtime, failcens) ~ treatment

res.strata <- get.strata.data(
  data.list       = data.list,
  ps_fml_covs     = ps.formula,
  v_arm           = NULL, # treat as a single-arm study since we want to borrow on both arms
  ctl_arm_level   = NULL,
  borrow_ctl_only = FALSE,
  nstrata         = nStrata,
  total_borrow    = nBorrow
)
data.list.PSIPP   <- res.strata$data.list
strata.list       <- res.strata$strata.list
res.strata$res.psrwe.borrow

fit.psipp <- aft.psipp(
  formula = fmla.psipp, 
  data.list = data.list.PSIPP
  , strata.list = strata.list
  , a0.strata = res.strata$a0.strata
  , dist = "weibull"
  , beta.mean = NULL, beta.sd = NULL
  , scale.mean = NULL, scale.sd = NULL
  , get.loglik = T,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
K <- attr(fit.psipp, 'data')$K
p <- attr(fit.psipp, 'data')$p
res.psipp <- suppressWarnings(
  fit.psipp %>%
    dplyr::select(c( paste0('treatment_stratum_', 1:K ),
                     paste0("scale_stratum_", 1:K))) %>%
    summarize_draws(mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
)
res.psipp

# compute log marginal likelihood
logml.psipp <- aft.logml.psipp(
  post.samples = fit.psipp
  , bridge.args = NULL
  , iter_warmup = warmup,
  iter_sampling = iter - warmup,
  chains = chains
)
summary(logml.psipp$bs)
summary(logml.psipp$bs.hist)
logml.psipp.val <- logml.psipp$logml

fit.list   <- list("Reference" = fit.ref, "PP" = fit.pp, "BHM" = fit.bhm,
                   "LEAP" = fit.leap, "PSIPP" = fit.psipp)
logml      <- c(logml.ref.val, logml.pp.val, logml.bhm.val, logml.leap.val, logml.psipp.val)

# estimated weights from each model-averaging method
loo.list  <- lapply(fit.list, function(f){
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
pbma_wts     <- as.numeric(pseudobma_weights(lpd_point, BB = FALSE))
pbma_BB_wts  <- as.numeric(pseudobma_weights(lpd_point, BB = TRUE)) # default is BB=TRUE
stacking_wts <- as.numeric(stacking_weights(lpd_point))
df.wts       <- data.frame(prior = names(fit.list),
                           "logml" = logml,
                           "BMA" = as.numeric(bridgesampling::post_prob(logml)),
                           "Pseudo-BMA" = pbma_wts,
                           "Pseudo-BMA w/ BB" = pbma_BB_wts,
                           "Stacking" = stacking_wts)
colnames(df.wts) <- c("prior", "logml", "BMA", "Pseudo-BMA", "Pseudo-BMA w/ BB", "Stacking")
df.wts %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))


# obtain posterior samples of predicted difference in 2-year survival probabilities
param_true_list <- readRDS(file = "~/Documents/UNC/Dissertation/super prior/Code/Sims_aft/R/trteff_true.rds")

# from each individual prior
surv.diff.2yr.list <- lapply(1:length(fit.list), function(i){
  prior <- names(fit.list)[i]
  f     <- fit.list[[i]]
  if( prior == "PSIPP" ){
    return(
      get.surv.diff.2yr.psipp(
        post.samples = f,
        trt.name = "treatment"
      )
    )
  }else{
    return(
      get.surv.diff.2yr(
        post.samples = f,
        trt.name = "treatment"
      )
    )
  }
})
df.estim <- sapply(1:length(surv.diff.2yr.list), function(i){
  s <- surv.diff.2yr.list[[i]]
  c(mean = mean(s), sd = sd(s), posterior::quantile2(s, probs = c(0.025, 0.975)))
})
df.estim        <- t(df.estim) %>% as.data.frame()
df.estim$method <- names(fit.list)
df.estim        <- df.estim[, c("method", colnames(df.estim)[which(colnames(df.estim) != "method")])]

# from each model-averaged prior
surv.diff.2yr.mtx      <- do.call(cbind, surv.diff.2yr.list)
surv.diff.2yr.list.avg <- lapply(colnames(df.wts)[-(1:2)], function(wts_name){
  wts     <- as.numeric( df.wts[, wts_name] )
  samples <- sample.model.avg.prior(wts = wts, samples.mtx = surv.diff.2yr.mtx)
  return(samples)
})
df.estim2 <- sapply(1:length(surv.diff.2yr.list.avg), function(i){
  s <- surv.diff.2yr.list.avg[[i]]
  c(mean = mean(s), sd = sd(s), posterior::quantile2(s, probs = c(0.025, 0.975)))
})
df.estim2        <- t(df.estim2) %>% as.data.frame()
df.estim2$method <- c("BMA", "Pseudo-BMA", "Pseudo-BMA w/ BB", "Stacking")
df.estim2        <- df.estim2[, c("method", colnames(df.estim2)[which(colnames(df.estim2) != "method")])]
df.estim         <- rbind(df.estim, df.estim2) %>% as.data.frame()

df.estim$param_true <- param_true_list$trteff.true
df.estim            <- df.estim %>%
  mutate(
    diff = mean - param_true,
    diff_sq = diff^2,
    ci.ind  = ifelse(param_true >= `q2.5` & param_true <= `q97.5`, 1, 0),
    log_var = ifelse(sd > 0, log(sd^2), NA)
  )
df.estim %>%
  mutate(across(where(is.numeric), \(x) round(x, 5)))

df.estim$method[order(df.estim$diff_sq)]
