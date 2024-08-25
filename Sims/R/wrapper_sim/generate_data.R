library(survival)
library(tidyverse)
library(psrwe)
library(dplyr)

mle.list <- readRDS(file = "Sims/R/mle_X_curr_hist.rds")
mle.curr <- mle.list$mle.curr
mle.hist <- mle.list$mle.hist
design.mtx.curr <- mle.list$design.mtx.curr
design.mtx.hist <- mle.list$design.mtx.hist
rm(mle.list)

## tdens: density function of survival times
get.cens.rate <- function(
    tdens, prop.cens = 0.30
) {
  optfun <- function(x) {
    f       <- function(t) { exp( log( tdens(t) ) - x * t ) }
    cenprop <- integrate(f, 0, Inf)$value
    abs(cenprop - prop.cens)
  }
  optimize(optfun, interval = c(0.01, 1000))$minimum
}

#' function to simulate time-to-event data from a Weibull AFT model
#' assume censoring time ~ exponential distribution and uniform accrual
#'
#' @param X             design matrix w/ column names being intercept/variable names
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current/historical data
#' @param accrual.max   maximum accrual time
#' For E1690 study, the accrual time was from February 1991 to June 1, 1995 (52/12 years),
#' and the presently analyzed database was current to September 1998
#' @param min.follow.up minimum follow up time for each subject, default is 3
#' @param censor.rate   censoring rate, default is log(0.95)/(-5)
#' @return data.frame with observed event time, event indicator, and covariates
get.weibull.surv <- function(
    X,
    mle,
    accrual.max = 52/12,
    min.follow.up = 3,
    censor.rate = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  X        <- as.matrix(X)
  n        <- nrow(X)
  sigma    <- exp(mle[1])  # MLE of scale parameter
  beta     <- mle[-1] # MLE of regression coefficients
  eta      <- X %*% beta

  # generate accrual time for each subject
  accrual.time     <- runif(n, min = 0, max = accrual.max)
  # administered censor time (the last subject enrolled in study will be observed for min.follow.up years)
  elapsed.time.max <- max(accrual.time) + min.follow.up
  # generate censoring time using exponential distribution with rate = censor.rate
  censor.time      <- rexp(n, rate = censor.rate)

  # generate failure/event time t from a Weibull distribution
  t <- rweibull(n, shape = 1/sigma, scale = exp(eta))

  # complete data observation time and event indicator
  observed.time <- pmin(t, censor.time)
  eventind      <- as.numeric(t <= censor.time)

  ## the complete data elapsed time
  elapsed.time  <- accrual.time + observed.time

  ## if the elapsed time for a subject is larger than the administered censor time,
  ## then the observed time is set to be (administered censor time - accrual time) and the event indicator is set to 0
  idx                <- (elapsed.time > elapsed.time.max)
  observed.time[idx] <- elapsed.time.max - (accrual.time[idx])
  eventind[idx]      <- 0

  if( "(Intercept)" %in% colnames(X) ) {
    X <- X[, !colnames(X) %in% "(Intercept)", drop = F]
  }
  df           <- cbind(observed.time, eventind, X)
  colnames(df) <- c('failtime', 'failcens', colnames(X))
  df           <- as.data.frame(df)
  return(df)
}


#' function to simulate data under the power prior (PP) assumption
#' generate current and historical data sets using mle.curr and bootstrap samples from X
#'
#' @param n             sample size of generated current data
#' @param n0            sample size of generated historical data
#' @param X             current data (E1690) design matrix, w/ column names being intercept/variable names
#' @param X0            historical data (E1684) design matrix, w/ column names being intercept/variable names
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current data
sim.PP <- function(
    n,
    n0,
    X             = design.mtx.curr,
    mle           = mle.curr,
    accrual.max   = 52/12,
    min.follow.up = 3,
    censor.rate   = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  # take bootstrap samples from current data
  idx       <- sample(1:nrow(X), size = n + n0, replace = T)
  X.all     <- X[idx, ] # combined design matrix for both generated current and generated historical data sets
  df.all    <- get.weibull.surv(
    X = X.all,
    mle = mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )
  data      <- df.all[(1:n), ]
  histdata  <- df.all[-(1:n), ]
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}


#' function to simulate data under the BHM assumption
#' generate current data set using mle.curr and bootstrap samples from X
#' generate historical data set using theta0 ~ N(mle.curr, mle.sd) and bootstrap samples from X
#'
#' @param n             sample size of generated current data
#' @param n0            sample size of generated historical data
#' @param X             current data (E1690) design matrix, w/ column names being intercept/variable names
#' @param X0            historical data (E1684) design matrix, w/ column names being intercept/variable names
#' @param mle.sd        standard deviation of the normal distribution from which the MLEs used to
#'                      generate the historical data are drawn
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current data
sim.BHM <- function(
    n,
    n0,
    X             = design.mtx.curr,
    mle.sd        = 0.25,
    mle           = mle.curr,
    accrual.max   = 52/12,
    min.follow.up = 3,
    censor.rate   = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  # take bootstrap samples from current data
  idx       <- sample(1:nrow(X), size = n + n0, replace = T)
  X.all     <- X[idx, ] # combined design matrix for both generated current and generated historical data sets
  data      <- get.weibull.surv(
    X = X.all[1:n, ],
    mle = mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )

  theta0    <- rnorm(n = length(mle), mean = mle, sd = mle.sd)
  histdata  <- get.weibull.surv(
    X = X.all[-(1:n), ],
    mle = theta0,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}


#' function to simulate data under the LEAP assumption
#' generate current data set using mle.curr
#' generate each observation in historical data set using mle.curr & X if tau = 1 and
#' using mle.hist & X0 if tau = 0, where tau ~ Bernoulli with probability = exch
#' (i.e., about exch*100 percent of subjects historical data are exchangeable with those in current data)
#'
#' @param n             sample size of generated current data
#' @param n0            sample size of generated historical data
#' @param X             current data (E1690) design matrix, w/ column names being intercept/variable names
#' @param X0            historical data (E1684) design matrix, w/ column names being intercept/variable names
#' @param exch          probability of subjects in historical data being exchangeable with those in
#'                      current data
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current data
#' @param mle0          MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on historical data
sim.LEAP <- function(
    n,
    n0,
    X             = design.mtx.curr,
    X0            = design.mtx.hist,
    exch          = 0.5,
    mle           = mle.curr,
    mle0          = mle.hist,
    accrual.max   = 52/12,
    min.follow.up = 3,
    censor.rate   = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  # take bootstrap samples from current data
  idx       <- sample(1:nrow(X), size = n, replace = T)
  data      <- get.weibull.surv(
    X = X[idx, ],
    mle = mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )

  # generate n0 Bernoulli random variables
  tau       <- rbinom(n0, size = 1, prob = exch)
  # generate historical data using mle & X for subjects with tau = 1
  idx01     <- sample(1:nrow(X), size = sum(tau), replace = T)
  histdata1 <- get.weibull.surv(
    X = X[idx01, ],
    mle = mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )
  # generate historical data using mle0 & X0 for subjects with tau = 0
  idx00     <- sample(1:nrow(X0), size = n0 - sum(tau), replace = T)
  histdata0 <- get.weibull.surv(
    X = X0[idx00, ],
    mle = mle0,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )
  histdata  <- rbind(histdata1, histdata0) %>%
    as.data.frame()
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}


#' function to simulate data under the unexchangeable assumption
#' generate current data set using mle.curr and bootstrap samples from X
#' generate historical data set using -mle.curr and bootstrap samples from X0
#'
#' @param n             sample size of generated current data
#' @param n0            sample size of generated historical data
#' @param X             current data (E1690) design matrix, w/ column names being intercept/variable names
#' @param X0            historical data (E1684) design matrix, w/ column names being intercept/variable names
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current data
sim.UNEXCH <- function(
    n,
    n0,
    X             = design.mtx.curr,
    X0            = design.mtx.hist,
    mle           = mle.curr,
    accrual.max   = 52/12,
    min.follow.up = 3,
    censor.rate   = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  # take bootstrap samples from current data
  idx       <- sample(1:nrow(X), size = n, replace = T)
  data      <- get.weibull.surv(
    X = X[idx, ],
    mle = mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )

  # take bootstrap samples from historical data
  idx0      <- sample(1:nrow(X0), size = n0, replace = T)
  histdata  <- get.weibull.surv(
    X = X0[idx0, ],
    mle = -mle,
    accrual.max = accrual.max,
    min.follow.up = min.follow.up,
    censor.rate = censor.rate
  )
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}


#' function to simulate data under the PSIPP assumption
#' generate current data set using mle.curr
#' generate historical data set using mle.curr and X if tau = 1 and using mle.hist and X0 if tau = 0,
#' where tau ~ Bernoulli with probability = exch
#'
#' @param n             sample size of generated current data
#' @param n0            sample size of generated historical data
#' @param X             current data (E1690) design matrix, w/ column names being intercept/variable names
#' @param X0            historical data (E1684) design matrix, w/ column names being intercept/variable names
#' @param mle           MLEs of scale parameter and regression coefficients from fitting Weibull AFT
#'                      model on current data
#' @param dev.mle.beta  deviation of the true regression coefficient values from `mle` for each stratum
#'
sim.PSIPP <- function(
    n,
    n0,
    X             = design.mtx.curr,
    X0            = design.mtx.hist,
    ps.formula    = ~ sex + cage + node_bin,
    nStrata       = 4,
    mle           = mle.curr,
    dev.mle.beta  = c(-0.2, -0.1, 0.1, 0.2),
    accrual.max   = 52/12,
    min.follow.up = 3,
    censor.rate   = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  # take bootstrap samples from current data
  idx       <- sample(1:nrow(X), size = n, replace = T)
  X.new     <- X[idx, ]
  # take bootstrap samples from historical data
  idx0      <- sample(1:nrow(X0), size = n0, replace = T)
  X0.new    <- X0[idx0, ]
  X.PS      <- rbind(X.new, X0.new)
  study     <- rep(c(1, 0), times = c(n, n0))
  df.PS     <- cbind(study = study, X.PS) %>%
    as.data.frame()

  # compute propensity scores (PS) via logistic regression and
  # obtain strata assignments using psrwe::psrwe_est() function
  res.PS.strata <- psrwe::psrwe_est(
    data = df.PS,
    ps_fml = as.formula( paste("study",  ps.formula) ),
    ps_method = "logistic",
    v_grp = "study",
    cur_grp_level = 1,
    v_arm = "treatment",
    ctl_arm_level = 0,
    stra_ctl_only = FALSE,
    nstrata = nStrata,
    trim_ab = "none"
  )
  df.PS.strata         <- res.PS.strata$data
  df.PS.strata$stratum <- as.integer( df.PS.strata$`_strata_` )

  # generate outcome within each stratum
  df.stratum.list <- lapply(1:nStrata, function(k){
    mle.stratum     <- mle
    mle.stratum[-1] <- mle.stratum[-1] + dev.mle.beta[k]
    X.stratum       <- df.PS.strata %>%
      filter(stratum == k)
    studyID.stratum <- X.stratum$study

    df.stratum      <- get.weibull.surv(
      X = X.stratum[, colnames(X)],
      mle = mle.stratum,
      accrual.max = accrual.max,
      min.follow.up = min.follow.up,
      censor.rate = censor.rate
    )
    df.stratum$study <- studyID.stratum
    return(df.stratum)
  })
  df.stratum <- do.call(rbind, df.stratum.list)

  data     <- df.stratum %>%
    filter(study == 1) %>%
    dplyr::select(-study)
  histdata <- df.stratum %>%
    filter(study == 0) %>%
    dplyr::select(-study)
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}


## source wrappers
wrapper.dir <- 'Analysis/R/wrapper_glm'
source(file.path(wrapper.dir, 'wrapper_surv.R'))
source(file.path(wrapper.dir, 'get_strata_data.R'))

#' function to transform data.list generated from the above sim.xxx functions (time-to-event data) into
#' counting process formats which can be used to fit different priors
#'
preprocess.data.list <- function(
    data.list,
    is.psipp   = TRUE,
    nbreaks    = 5,
    ps.formula = ~ sex + cage + node_bin,
    nStrata    = 4
){
  probs   <- 1:nbreaks / nbreaks
  breaks  <- data.list[[1]] %>%
    filter(failcens == 1) %>%
    reframe(quant = quantile(failtime, probs = probs)) %>%
    unlist
  breaks  <- breaks[-nbreaks]

  if( is.psipp ){
    nBorrow    <- nrow(data.list[[2]])
    res.strata <- get.strata.data(
      data.list       = data.list,
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

    ## force historical data to have the same number of levels for the interval variable as the current data
    if( any( !levels(curr.poisson$interval) %in% levels(hist.poisson$interval) ) ){
      hist.poisson$interval <- factor(hist.poisson$interval, levels = levels(curr.poisson$interval))
    }

    data.list      <- list(currdata = curr.poisson, histdata = hist.poisson)
    offset.list    <- list(offs = curr.poisson$log_exposue, offs0 = hist.poisson$log_exposue)
    strata.list    <- lapply(data.list, function(l){
      l[, 'ps_strata']
    })
    res <- list(
      data.list   = data.list,
      offset.list = offset.list,
      strata.list = strata.list,
      a0.strata   = res.strata$a0.strata,
      res.strata  = res.strata,
      breaks      = breaks
    )

  }else{
    ## transform data into counting process format
    covariate.name <- c("treatment", "sex", "cage", "node_bin")
    curr.poisson   <- get.counting.process.data(
      data = data.list[[1]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
      covariate.name = covariate.name
    )
    hist.poisson   <- get.counting.process.data(
      data = data.list[[2]], breaks = breaks, failtime.name = "failtime", failind.name = "failcens",
      covariate.name = covariate.name
    )

    ## force historical data to have the same number of levels for the interval variable as the current data
    if( any( !levels(curr.poisson$interval) %in% levels(hist.poisson$interval) ) ){
      hist.poisson$interval <- factor(hist.poisson$interval, levels = levels(curr.poisson$interval))
    }

    data.list      <- list(currdata = curr.poisson, histdata = hist.poisson)
    offset.list    <- list(offs = curr.poisson$log_exposue, offs0 = hist.poisson$log_exposue)
    res <- list(
      data.list   = data.list,
      offset.list = offset.list,
      breaks      = breaks
    )

  }
  return(res)
}
