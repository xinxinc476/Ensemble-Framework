# check the poisson likelihood computation is equivalent (or differs only by a constant) to
# the likelihood of fitting the piece-wise exponential model (PWE model)
library(tidyverse)
library(survival)
library(hdbayes)
library(posterior)

# load data
curr <- E1690

# replace 0 failure times with 0.50 days
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 10 subjects w/ failtime = 0

# Center and scale age based on current data
age_stats <- with(curr, c('mean' = mean(age), 'sd' = sd(age)))
curr$cage <- ( curr$age - age_stats['mean'] ) / age_stats['sd']

# assume 5 intervals
# obtain cut points for intervals
nbreaks <- 5
probs   <- 1:nbreaks / nbreaks
breaks  <- curr %>%
  filter(failcens == 1) %>%
  reframe(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks  <- breaks[-nbreaks]

# compute the poisson log-likelihood
family <- poisson("log")
fmla   <- failcens ~ interval + treatment + sex + cage + node_bin

wrapper.dir <- 'Analysis/R/wrapper_glm'
source(file.path(wrapper.dir, 'wrapper_surv.R'))

# transform data into counting process format
covariate.name <- c("treatment", "sex", "cage", "node_bin")
curr.poisson   <- get.counting.process.data(
  data = curr, breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
  covariate.name = covariate.name
)

data.list      <- list(currdata = curr.poisson)
offset.list    <- list(offs = curr.poisson$log_exposue)

# fit a non-informative reference prior
d <- glm.post(
  formula = fmla, family = family, data.list = data.list,
  offset.list = offset.list,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

post.samples <- rbind(colMeans(d), d)
loglik       <- get.loglik.glm(post.samples)
loglik.draws <- rowSums( loglik ) # log-likelihood from poisson model

# function to compute log-likelihood under the PWE model
pwe.loglik <- function(beta, log_hazard) {
  y      <- curr$failtime
  X      <- model.matrix(~treatment + sex + cage + node_bin, curr)
  # remove intercept term
  X      <- X[, -1]
  fail   <- curr$failcens

  # compute interval index for which each observed time y belongs to
  intindx <- sapply(y, function(t){
    findInterval(t, c(0, breaks, Inf), left.open = TRUE)
  })

  # compute baseline hazards (lambda)
  log_hazard <- log_hazard[1] + c(0, log_hazard[2:nbreaks])
  lambda     <- as.numeric(exp(log_hazard))
  # compute linear predictor
  eta        <- as.numeric( (X %*% beta)[, 1] )

  # get hazard corresponding to failure / censoring time
  lambda_j   <- lambda[intindx]

  # compute cumulative baseline hazard at each interval
  breaks.new  <- c(0, breaks, Inf)
  cumblhaz    <- cumsum( lambda[1:(nbreaks-1)] * ( breaks.new[2:nbreaks] - breaks.new[1:(nbreaks-1)] ) )
  cumblhaz    <- c(0, cumblhaz)

  # compute cumulative hazard for each observation
  cumhaz      <- lambda_j * (y - breaks.new[intindx]) + cumblhaz[intindx]
  cumhaz      <- cumhaz * exp(eta)

  # Log likelihood = sum( event_ind * log(hazard) - cumhaz ),
  # where log(hazard) = log( lambda * exp(eta) ) = log(lambda) + eta
  loglik <- sum( fail * ( log(lambda_j) + eta ) - cumhaz )
  return(loglik)
}

loglik.draws.pwe <- sapply(1:nrow(post.samples), function(i){
  param <- post.samples[i, ]
  beta.i <- suppressWarnings(
    as.numeric(param[covariate.name])
  )
  log_hazard.i <- suppressWarnings(
    as.numeric(param[2:6])
  )
  return(
    pwe.loglik(beta = beta.i ,
               log_hazard = log_hazard.i)
  )
}) # log-likelihood from PWE model

diff.loglik <- loglik.draws - loglik.draws.pwe
table(diff.loglik) # differ by a constant

all.equal(
  diff.loglik[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE
