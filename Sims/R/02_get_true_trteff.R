## compute the true value of the estimand: difference in 2-year relapse-free survival (RFS) probability,
## for treated v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)
library(tidyverse)
library(doParallel)
library(MCMCpack)

mle.list <- readRDS(file = "Sims/R/mle_X_curr_hist.rds")
mle.curr <- mle.list$mle.curr
mle.hist <- mle.list$mle.hist
design.mtx.curr <- mle.list$design.mtx.curr
design.mtx.hist <- mle.list$design.mtx.hist
rm(mle.list)

# wrapper function to compute the true value of the estimand
get.true.trteff <- function(nevents = 1000){
  X        = as.matrix(design.mtx.curr)
  theta    = mle.curr
  n        = nevents
  # take bootstrap samples from X
  idx      = sample(1:nrow(X), size = n, replace = T)
  X.new    = X[idx, ]

  # generate treatment indicator
  X.trt                = X.new
  X.trt[, "treatment"] = 1 ## assume all subjects are in the treated group
  X.ctl                = X.new
  X.ctl[, "treatment"] = 0 ## assume all subjects are in the control group
  
  sigma    = as.numeric( theta[1] )  # MLE of sigma
  beta     = as.numeric( theta[-1] ) # MLE of regression coefficients
  eta.trt  = as.numeric( X.trt %*% beta )
  eta.ctl  = as.numeric( X.ctl %*% beta )
  
  # compute 2-year RFS probability for each subject,
  # assuming the RFS time follows a Weibull(alpha, lambda) distribution
  # shape = alpha = 1/sigma, scale = exp(-lambda / alpha), where lambda = -x'beta/sigma
  surv.prob.trt = pweibull(2, shape = 1/sigma, scale = exp(eta.trt), lower.tail = FALSE)
  surv.prob.ctl = pweibull(2, shape = 1/sigma, scale = exp(eta.ctl), lower.tail = FALSE)

  ## use Bayesian bootstrap method
  ## sample from Dirichlet(1, 1, .., 1) distribution
  omega.trt  = as.numeric( MCMCpack::rdirichlet(n = 1, alpha = rep(1, nrow(X.trt))) )
  omega.ctl  = as.numeric( MCMCpack::rdirichlet(n = 1, alpha = rep(1, nrow(X.ctl))) )
  surv.trt   = sum( omega.trt * surv.prob.trt )
  surv.ctl   = sum( omega.ctl * surv.prob.ctl)
  return(
    surv.trt - surv.ctl
  )
}

set.seed(100)
trteff.true <- get.true.trteff(nevents = 1000000)

saveRDS(
  list(
    trteff.true = trteff.true
  ), file = "Sims/R/trteff_true.rds"
)
