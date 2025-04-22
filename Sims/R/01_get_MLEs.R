# obtain MLEs from Weibull AFT models using current and historical data
library(survival)
library(hdbayes)
library(tidyverse)

## load data
hist <- E1684
curr <- E1690

## replace 0 failure times with 0.50 days
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 1 subject w/ failtime = 0
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 10 subjects w/ failtime = 0

## Center and scale age based on current data
age_stats <- with(curr, c('mean' = mean(age), 'sd' = sd(age)))
hist$cage <- ( hist$age - age_stats['mean'] ) / age_stats['sd']
curr$cage <- ( curr$age - age_stats['mean'] ) / age_stats['sd']

# fit Weibull models to current and historical data
fmla     <- Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
fit.curr <- survreg(formula = fmla, data = curr, dist = "weibull")
fit.hist <- survreg(formula = fmla, data = hist, dist = "weibull")

# estimated log(scale) and regression coefficients from Weibull AFT models
mle.curr <- c(scale = fit.curr$scale, fit.curr$coefficients)
mle.hist <- c(scale = fit.hist$scale, fit.hist$coefficients)

# design matrix from current and historical data sets
X        <- model.matrix(fmla, curr)
X0       <- model.matrix(fmla, hist)

# save MLEs
saveRDS(
  list(
    mle.curr = mle.curr,
    mle.hist = mle.hist,
    design.mtx.curr = X,
    design.mtx.hist = X0
  ),
  file = "Sims/R/mle_X_curr_hist.rds"
)
