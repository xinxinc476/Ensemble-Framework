library(survival)
library(dplyr)

source("wrapper/wrapper_sim.R")

#' function to simulate data under the unexchangeable case
#'
#' @param prop.cens     desired censoring proportion
#' @param nevents       desired number of events in simulated current data
#' @param n0events      desired number of events in simulated historical data
#' @param X             design matrix w/ column names being "(intercept)", "treatment", and covariate
#'                      names (if any). Bootstrap samples from the covariates of `X` will be used as
#'                      the covariates for simulated current data. Defaults to design matrix from
#'                      E1690 data.
#' @param X0            design matrix w/ column names being "(intercept)", "treatment", and covariate
#'                      names (if any). Bootstrap samples from the covariates of `X0` will be used as
#'                      the covariates for simulated historical data. Defaults to design matrix from
#'                      E1684 data.
#' @param theta         a vector of scale parameter and regression coefficients, e.g., the MLEs
#'                      from fitting a Weibull AFT model on current/historical data. `theta` will be
#'                      used for simulating current data.
#' @param theta0        a vector of scale parameter and regression coefficients, e.g., the MLEs
#'                      from fitting a Weibull AFT model on current/historical data. `theta0` will be
#'                      used for simulating historical data.
#' @param trt.prob      probability of being assigned to the treated group. The treatment indicator
#'                      for simulated data will be generated using Bernoulli distribution with probability
#'                      being `trt.prob`. Defaults to 0.5.
#'
sim.UNEXCH <- function(
    prop.cens,
    nevents,
    n0events,
    X             = design.mtx.curr,
    X0            = design.mtx.hist,
    theta         = mle.curr,
    theta0        = 2 * mle.hist,
    trt.prob      = 0.5
){
  # equivalent to simulate under LEAP assumption with `exch` = 0
  data      <- get.weibull.surv(
    X = X,
    theta = theta,
    prop.cens = prop.cens,
    nevents = nevents,
    trt.prob = trt.prob,
    bootstrap = TRUE
  )

  # generate historical data using theta0 & X0
  histdata  <- get.weibull.surv(
    X = X0,
    theta = theta0,
    prop.cens = prop.cens,
    nevents = n0events,
    trt.prob = trt.prob,
    bootstrap = TRUE
  )

  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}
