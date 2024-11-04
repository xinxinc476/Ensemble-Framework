# check the poisson likelihood computation is equivalent (or differs only by a
# constant) to the likelihood of the piece-wise exponential model (PWE model)
library(tidyverse)
library(survival)
library(hdbayes)
library(posterior)

# load data
curr <- E1690
hist <- E1684

# replace 0 failure times with 0.50 days
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) 
# 10 subjects w/ failtime = 0
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) 
# 1 subject w/ failtime = 0

# Center and scale age based on current data
age_stats <- with(curr, c('mean' = mean(age), 'sd' = sd(age)))
curr$cage <- ( curr$age - age_stats['mean'] ) / age_stats['sd']
hist$cage <- ( hist$age - age_stats['mean'] ) / age_stats['sd']

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

################################# check for non-PSIPP prior ###############################

wrapper.dir <- "R/wrapper_glm"
source(file.path(wrapper.dir, "wrapper_surv.R"))

fmla   <- failcens ~ interval + treatment + sex + cage + node_bin

# transform data into counting process format
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

# function to compute log-likelihood matrix under the PWE model (for non-PSIPP prior)
pwe.loglik <- function(beta, log_hazard, data = curr) {
  y      <- data$failtime
  X      <- model.matrix(~treatment + sex + cage + node_bin, data)
  # remove intercept term
  X      <- X[, -1]
  fail   <- data$failcens
  
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
  
  # Log likelihood = event_ind * log(hazard) - cumhaz,
  # where log(hazard) = log( lambda * exp(eta) ) = log(lambda) + eta
  loglik <- fail * ( log(lambda_j) + eta ) - cumhaz 
  return(loglik)
}


################################# reference prior ###############################

# fit a non-informative reference prior
d <- glm.post(
  formula = fmla, family = family, data.list = data.list,
  offset.list = offset.list,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

loglik       <- get.loglik.glm(d) # log-likelihood matrix from poisson model
loglik.draws <- rowSums( loglik ) # sum of log-likelihood from poisson model

loglik.pwe <- sapply(1:nrow(d), function(i){
  param <- d[i, ]
  beta.i <- suppressWarnings(
    as.numeric(param[covariate.name])
  )
  log_hazard.i <- suppressWarnings(
    as.numeric(param[2:(nbreaks+1)])
  )
  return(
    pwe.loglik(beta = beta.i ,
               log_hazard = log_hazard.i)
  )
}) 
loglik.pwe       <- t(loglik.pwe) # log-likelihood from PWE model
loglik.pwe.draws <- rowSums( loglik.pwe ) # sum of log-likelihood from PWE model

diff.loglik <- loglik.draws - loglik.pwe.draws # difference in sum of log-likelihood for each MCMC iteration
table(round(diff.loglik, 9)) # differ by a constant up to 9 decimal places

all.equal(
  diff.loglik[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE

loglik2 <- loglik %r-% (curr.poisson$log_exposue * curr.poisson$failcens)
all.equal(rowSums( loglik2 ), loglik.pwe.draws ) # TRUE


################################# power prior (PP) ###############################

# fit PP w/ a0 = 0.5
d.pp <- glm.pp(
  formula = fmla, family = family, data.list = data.list,
  offset.list = offset.list,
  a0 = 0.5,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

loglik.pp       <- get.loglik.glm(d.pp) # log-likelihood matrix from poisson model
loglik.pp.draws <- rowSums( loglik.pp ) # sum of log-likelihood from poisson model

loglik.pp.pwe <- sapply(1:nrow(d.pp), function(i){
  param <- d.pp[i, ]
  beta.i <- suppressWarnings(
    as.numeric(param[covariate.name])
  )
  log_hazard.i <- suppressWarnings(
    as.numeric(param[2:(nbreaks+1)])
  )
  return(
    pwe.loglik(beta = beta.i ,
               log_hazard = log_hazard.i)
  )
}) 
loglik.pp.pwe       <- t(loglik.pp.pwe) # log-likelihood from PWE model
loglik.pp.pwe.draws <- rowSums( loglik.pp.pwe ) # sum of log-likelihood from PWE model

diff.loglik.pp <- loglik.pp.draws - loglik.pp.pwe.draws
table(round(diff.loglik.pp, 9)) # differ by a constant up to 9 decimal places

all.equal(
  diff.loglik.pp[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE

loglik2.pp <- loglik.pp %r-% (curr.poisson$log_exposue * curr.poisson$failcens)
all.equal(rowSums( loglik2.pp ), loglik.pp.pwe.draws ) # TRUE


################################### BHM ##################################

# fit BHM
d.bhm <- glm.bhm(
  formula = fmla, family = family, data.list = data.list,
  offset.list = offset.list,
  meta.sd.sd = 1,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

loglik.bhm       <- get.loglik.glm(d.bhm) # log-likelihood matrix from poisson model
loglik.bhm.draws <- rowSums( loglik.bhm ) # sum of log-likelihood from poisson model

loglik.bhm.pwe <- sapply(1:nrow(d.bhm), function(i){
  param <- d.bhm[i, ]
  beta.i <- suppressWarnings(
    as.numeric(param[covariate.name])
  )
  log_hazard.i <- suppressWarnings(
    as.numeric(param[2:(nbreaks+1)])
  )
  return(
    pwe.loglik(beta = beta.i ,
               log_hazard = log_hazard.i)
  )
}) 
loglik.bhm.pwe       <- t(loglik.bhm.pwe) # log-likelihood from PWE model
loglik.bhm.pwe.draws <- rowSums( loglik.bhm.pwe ) # sum of log-likelihood from PWE model

diff.loglik.bhm <- loglik.bhm.draws - loglik.bhm.pwe.draws
table(round(diff.loglik.bhm, 9)) # differ by a constant up to 9 decimal places

all.equal(
  diff.loglik.bhm[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE

loglik2.bhm <- loglik.bhm %r-% (curr.poisson$log_exposue * curr.poisson$failcens)
all.equal(rowSums( loglik2.bhm ), loglik.bhm.pwe.draws ) # TRUE


################################### LEAP #################################

# fit LEAP w/ K = 2
K           <- 2
offset.list.leap <- lapply(offset.list, function(l){
  matrix(rep(l, K), ncol = K)
})
d.leap <- glm.leap(
  formula = fmla, family = family, data.list = data.list,
  offset.list = offset.list.leap,
  K = K,
  gamma.upper = 1,
  gamma.lower = 0,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

loglik.leap       <- get.loglik.glm(d.leap) # log-likelihood matrix from poisson model
loglik.leap.draws <- rowSums( loglik.leap ) # sum of log-likelihood from poisson model

loglik.leap.pwe <- sapply(1:nrow(d.leap), function(i){
  param <- d.leap[i, ]
  beta.i <- suppressWarnings(
    as.numeric(param[covariate.name])
  )
  log_hazard.i <- suppressWarnings(
    as.numeric(param[2:(nbreaks+1)])
  )
  return(
    pwe.loglik(beta = beta.i ,
               log_hazard = log_hazard.i)
  )
}) 
loglik.leap.pwe       <- t(loglik.leap.pwe) # log-likelihood from PWE model
loglik.leap.pwe.draws <- rowSums( loglik.leap.pwe ) # sum of log-likelihood from PWE model

diff.loglik.leap <- loglik.leap.draws - loglik.leap.pwe.draws
table(round(diff.loglik.leap, 9)) # differ by a constant up to 9 decimal places

all.equal(
  diff.loglik.leap[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE

loglik2.leap <- loglik.leap %r-% (curr.poisson$log_exposue * curr.poisson$failcens)
all.equal(rowSums( loglik2.leap ), loglik.leap.pwe.draws ) # TRUE


##################################### check for PSIPP ##################################

wrapper.dir <- "R/wrapper_glm"
source(file.path(wrapper.dir, "get_strata_data.R"))
source(file.path(wrapper.dir, "glm_stratified_pp.R"))
source(file.path(wrapper.dir, "glm_logml_stratified_pp.R"))

fmla.psipp <- failcens ~ interval + treatment
ps.formula <- ~ sex + cage + node_bin
nStrata    <- 4 ## number of strata for analysis using PSIPP
nBorrow    <- nrow(hist)

# add subject ID
curr$ID <- 1:nrow(curr)
hist$ID <- 1:nrow(hist)

res.strata <- get.strata.data(
  data.list       = list(curr, hist),
  ps_fml_covs     = ps.formula,
  #v_arm           = "treatment",
  #ctl_arm_level   = 0,
  v_arm           = NULL,
  ctl_arm_level   = NULL,
  borrow_ctl_only = FALSE,
  nstrata         = nStrata,
  total_borrow    = nBorrow
)
data.list.psipp <- res.strata$data.list
strata.list     <- res.strata$strata.list
data.list.psipp <- lapply(1:length(data.list.psipp), function(i){
  df = data.list.psipp[[i]]
  df$ps_strata = strata.list[[i]]
  return(df)
})

# transform data into counting process format
covariate.name <- c("treatment", "sex", "cage", "node_bin", "ps_strata", "ID")
curr.poisson   <- get.counting.process.data(
  data = data.list.psipp[[1]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
  covariate.name = covariate.name
)
hist.poisson   <- get.counting.process.data(
  data = data.list.psipp[[2]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
  covariate.name = covariate.name
)
data.list.psipp.poisson   <- list(currdata = curr.poisson, histdata = hist.poisson)
offset.list.psipp.poisson <- list(offs = curr.poisson$log_exposue, offs0 = hist.poisson$log_exposue)
strata.list.poisson       <- lapply(data.list.psipp.poisson, function(l){
  l[, 'ps_strata']
})

# function to compute log-likelihood matrix under the PWE model (for PSIPP)
pwe.loglik.psipp <- function(beta, log_hazard, data = data.list.psipp[[1]]) {
  y      <- data$failtime
  X      <- data$treatment # treatment indicator
  fail   <- data$failcens
  
  # compute interval index for which each observed time y belongs to
  intindx <- sapply(y, function(t){
    findInterval(t, c(0, breaks, Inf), left.open = TRUE)
  })
  
  # compute log-likelihood matrix for each stratum
  loglik.psipp <- sapply(1:nStrata, function(k){
    y_stratum          <- y[data$ps_strata == k]
    fail_stratum       <- fail[data$ps_strata == k]
    X_stratum          <- X[data$ps_strata == k]
    intindx_stratum    <- intindx[data$ps_strata == k]
    
    # compute baseline hazards (lambda)
    log_hazard_stratum <- as.numeric(log_hazard[, paste0(interval.names, '_stratum_', k)])
    log_hazard_stratum <- log_hazard_stratum[1] + c(0, log_hazard_stratum[2:nbreaks])
    lambda_stratum     <- exp(log_hazard_stratum) ## baseline hazards
    
    # get hazard corresponding to failure / censoring time
    lambda_j        <- lambda_stratum[intindx_stratum]
    
    # compute cumulative baseline hazard at each interval
    breaks.new  <- c(0, breaks, Inf)
    cumblhaz    <- cumsum( lambda_stratum[1:(nbreaks-1)] * ( breaks.new[2:nbreaks] - breaks.new[1:(nbreaks-1)] ) )
    cumblhaz    <- c(0, cumblhaz)
    
    # compute linear predictor
    eta         <- X_stratum * as.numeric(beta[k])
    
    # compute cumulative hazard for each observation
    cumhaz      <- lambda_j * (y_stratum - breaks.new[intindx_stratum]) + cumblhaz[intindx_stratum]
    cumhaz      <- cumhaz * exp(eta)
    
    # Log likelihood = event_ind * log(hazard) - cumhaz,
    # where log(hazard) = log( lambda * exp(eta) ) = log(lambda) + eta
    return( fail_stratum * ( log(lambda_j) + eta ) - cumhaz )
  })
  loglik.psipp <- do.call(c, loglik.psipp)
  return(loglik.psipp)
}

d.psipp <- glm.stratified.pp(
  formula = fmla.psipp, family = family, 
  data.list = data.list.psipp.poisson,
  strata.list = strata.list.poisson,
  a0.strata = res.strata$a0.strata,
  offset.list = offset.list.psipp.poisson,
  is.sorted = FALSE,
  iter_warmup = 2000,
  iter_sampling = 20000,
  chains = 1
)

loglik.psipp       <- get.loglik.glm.stratified.pp(d.psipp) # log-likelihood matrix from poisson model
loglik.psipp.draws <- rowSums( loglik.psipp ) # sum of log-likelihood from poisson model

interval.names    <- paste0("interval", levels(curr.poisson$interval))
interval.names[1] <- c("(Intercept)")
loglik.psipp.pwe  <- sapply(1:nrow(d.psipp), function(i){
  param <- d.psipp[i, ]
  beta.i <- suppressWarnings(
    param[paste0( "treatment", '_stratum_', 1:nStrata )]
  )
  
  log_hazard.i <- suppressWarnings(
    param[paste0( interval.names, '_stratum_', rep(1:nStrata, each = nbreaks) )]
  )
  return(
    pwe.loglik.psipp(beta = beta.i, log_hazard = log_hazard.i)
  )
}) 
loglik.psipp.pwe       <- t(loglik.psipp.pwe) # log-likelihood from PWE model
loglik.psipp.pwe.draws <- rowSums( loglik.psipp.pwe ) # sum of log-likelihood from PWE model

diff.loglik.psipp <- loglik.psipp.draws - loglik.psipp.pwe.draws
table(round(diff.loglik.psipp, 9)) # differ by a constant up to 9 decimal places
all.equal(diff.loglik.psipp[1], diff.loglik[1]) # true

all.equal(
  diff.loglik[1],
  sum(curr.poisson$log_exposue * curr.poisson$failcens)
) # TRUE

idx <- order(data.list.psipp.poisson$currdata$ID)
loglik.psipp2 <- loglik.psipp %r-% (data.list.psipp.poisson$currdata$log_exposue[idx] * data.list.psipp.poisson$currdata$failcens[idx])
all.equal(rowSums( loglik.psipp2 ), loglik.psipp.pwe.draws ) # TRUE

