## compute the true value of the estimand: difference in survival (or relapse-free survival) probability at 2 years for treated
## v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)

## simulated current data follow the same distribution under the following cases: UNEXCH, LEAP, BHM, and PP

library(tidyverse)
library(doParallel)
library(MCMCpack)

mle.list <- readRDS(file = "Sims_aft/R/mle_X_curr_hist.rds")
mle.curr <- mle.list$mle.curr
mle.hist <- mle.list$mle.hist
design.mtx.curr <- mle.list$design.mtx.curr
design.mtx.hist <- mle.list$design.mtx.hist
rm(mle.list)

# wrapper function to obtain compute the true value of the estimand under UNEXCH, LEAP, BHM, or PP
get.true.trteff <- function(iter, nevents = 1000){
  X        = as.matrix(design.mtx.curr)
  theta    = mle.curr
  n        = nevents
  # take bootstrap samples from X
  idx      = sample(1:nrow(X), size = n, replace = T)
  X.new    = X[idx, ]
  # generate treatment indicator
  X.new[, which(colnames(X.new) == "treatment")] <- rbinom(n, 1, 0.5)
  
  sigma    = as.numeric( theta[1] )  # MLE of sigma
  beta     = as.numeric( theta[-1] ) # MLE of regression coefficients
  eta      = as.numeric( X.new %*% beta )
  
  # compute 2-year survival probability for each subject,
  # assuming the survival time follows a Weibull(alpha, lambda) distribution
  # shape = alpha = 1/sigma, scale = exp(-lambda / alpha), where lambda = -x'beta/sigma
  surv.prob  = pweibull(2, shape = 1/sigma, scale = exp(eta), lower.tail = FALSE)
  
  id.trt     = which( X.new[, "treatment"] == 1 )
  id.ctl     = which( X.new[, "treatment"] == 0 )
  
  ## use Bayesian bootstrap to compute 2-year survival probability for subjects in the treatment and control arms
  ## sample from dirichlet(1, 1, .., 1) distribution
  omega.trt  = as.numeric( MCMCpack::rdirichlet(n = 1, alpha = rep(1, length(id.trt))) )
  omega.ctl  = as.numeric( MCMCpack::rdirichlet(n = 1, alpha = rep(1, length(id.ctl))) )
  surv.trt   = sum( omega.trt * surv.prob[id.trt] )
  surv.ctl   = sum( omega.ctl * surv.prob[id.ctl])
  return(
    surv.trt - surv.ctl
  )
}

nsims <- 1000000
cl    <- makeCluster(10)
registerDoParallel(cl)
res   <- foreach(i = 1:nsims, .combine=c) %dopar%
  get.true.trteff(iter = i, nevents = 1000)
trteff.true <- mean(res)

saveRDS(
  list(
    trteff.true = trteff.true
    #, trteff.true.psipp = trteff.true.psipp
  ), file = "Sims_aft/R/trteff_true.rds"
)


## compute the true value of the estimand under PSIPP
theta      <- mle.curr[c("logscale", "(Intercept)", "treatment" )]
sigma      <- as.numeric( exp(theta[1]) )  ## MLE of sigma
nStrata    <- 3
drift.beta <- c(-0.1, 0, 0.1)
## predicted 2-year survival probability for untreated for each stratum
surv.ctl   <- lapply(1:nStrata, function(k){
  theta.stratum     = theta
  theta.stratum[-1] = theta.stratum[-1] + drift.beta[k]
  S_ctl_stratum     = pweibull(2, shape = 1/sigma,
                               scale = exp(theta.stratum["(Intercept)"]),
                               lower.tail = FALSE)
  return(S_ctl_stratum)
}
)
surv.ctl   <- unlist(surv.ctl)

## predicted 2-year survival probability for treated for each stratum
surv.trt   <- lapply(1:nStrata, function(k){
  theta.stratum     = theta
  theta.stratum[-1] = theta.stratum[-1] + drift.beta[k]
  S_trt_stratum     = pweibull(2, shape = 1/sigma,
                               scale = exp(theta.stratum["(Intercept)"] + theta.stratum["treatment"]),
                               lower.tail = FALSE)
  return(S_trt_stratum)
}
)
surv.trt   <- unlist(surv.trt)
trteff.true.psipp <- mean(surv.trt) - mean(surv.ctl)

saveRDS(list(trteff.true = trteff.true, trteff.true.psipp = trteff.true.psipp),
        file = "Sims/R/trteff_true.rds")
