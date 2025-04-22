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


#' function to compute the optimal exponential rate to achieve the desired censoring proportion
#' assume event time follows a Weibull distribution as specified in [stats::rweibull()]
#' assume censoring time follows a exponential distribution
#'
#' @param sigma         scale parameter from AFT model: Y = log(T) = mu + sigma * error term
#' @param mu            parameter from AFT model: Y = log(T) = mu + sigma * error term
#' @param prop.cens     desired censoring proportion
#'
#' @examples
#' sigma <- 1.5
#' mu <- 0.3
#' t <- rweibull(10000, shape = 1/sigma, scale = exp(mu))
#' x <- get.cens.rate(sigma = sigma, mu = mu, prop.cens = 0.35)
#' c <- rexp(10000, rate = x)
#' mean(c < t) # should be around 0.35
#'
get.cens.rate <- function(
    sigma,
    mu,
    prop.cens
) {
  optfun <- function(x) {
    f       <- function(t){
      stats::pweibull(t, shape = 1/sigma, scale = exp(mu), lower.tail = F, log.p = F) * stats::dexp(t, rate = x, log = F)
    }
    censprop <- integrate(f, 0, Inf)$value
    abs(censprop - prop.cens)
  }
  stats::optimize(optfun, interval = c(0.001, 1000))$minimum
}


#' function to simulate time-to-event data from a Weibull accelerated failure time (AFT) model
#' assume censoring time ~ exponential distribution with the rate parameter being determined based on
#' the desired censoring proportion.
#'
#'
#' @param X             design matrix w/ column names being "(intercept)", "treatment", and covariate
#'                      names (if any). Bootstrap samples from the covariates of `X` will be used as
#'                      the covariates for simulated data.
#' @param theta         a vector of scale parameter and regression coefficients, e.g., the MLEs
#'                      from fitting a Weibull AFT model on current/historical data.
#' @param prop.cens     desired censoring proportion
#' @param nevents       desired number of events
#' @param trt.prob      probability of being assigned to the treated group. The treatment indicator
#'                      for simulated data will be generated using Bernoulli distribution with probability
#'                      being `trt.prob`. Defaults to 0.5.
#' @param bootstrap     whether to take bootstrap samples from `X` with replacement and generate treatment
#'                      indicators. If FALSE, `X` will be used as the design matrix for new data, and
#'                      arguments like `nevents` and `trt.prob` are not used. Defaults to TRUE.
#'
#' @return
#'  a data.frame with observed event time, event indicator, treatment indicator, and covariates.
#'
#' @examples
#' df <- get.weibull.surv(
#'   X         = design.mtx.curr,
#'   theta     = mle.curr,
#'   prop.cens = 0.2,
#'   nevents   = 200,
#'   trt.prob  = 0.5,
#'   bootstrap = TRUE
#' )
#' mean(1 - df$failcens) # should be around 0.2
#' # plot KM curve
#' fit  <- survfit(Surv(failtime, failcens) ~ treatment , data = df)
#' survminer::ggsurvplot(fit, data = df)
#'
get.weibull.surv <- function(
    X,
    theta,
    prop.cens,
    nevents,
    trt.prob = 0.5,
    bootstrap = TRUE
){
  X        <- as.matrix(X)
  if( bootstrap ){
    n        <- ceiling( nevents / (1 - prop.cens) )
    # take bootstrap samples from X
    idx      <- sample(1:nrow(X), size = n, replace = T)
    X.new    <- X[idx, ]

    # generate treatment indicator
    X.new[, which(colnames(X.new) == "treatment")] <- rbinom(n, 1, trt.prob)
  }else{
    n     <- nrow(X)
    X.new <- X
  }

  sigma    <- as.numeric( theta[1] )  # MLE of sigma
  beta     <- as.numeric( theta[-1] ) # MLE of regression coefficients
  eta      <- as.numeric( X.new %*% beta )

  # generate event time t from a Weibull distribution
  t        <- rweibull(n, shape = 1/sigma, scale = exp(eta))

  if ( prop.cens > 0 ){
    # compute the appropriate censoring rate to achieve the desired censoring proportion
    censor.rate <- sapply(eta, function(z){
      get.cens.rate(sigma = sigma, mu = z, prop.cens = prop.cens)
    })
    # generate censoring time using exponential distribution with rate = censor.rate
    censor.time <- rexp(n = n, rate = censor.rate)

    # compute observation time and event indicator
    observed.time <- pmin(t, censor.time)
    eventind      <- as.numeric(t <= censor.time)
  }else{
    observed.time <- t
    eventind      <- rep(1, n)
  }

  if( "(Intercept)" %in% colnames(X.new) ) {
    X.new <- X.new[, !colnames(X.new) %in% "(Intercept)", drop = F]
  }
  df           <- cbind(observed.time, eventind, X.new)
  colnames(df) <- c('failtime', 'failcens', colnames(X.new))
  df           <- as.data.frame(df)
  return(df)
}
