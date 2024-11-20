library(survival)
library(dplyr)

source("R/wrapper_aft/wrapper_sim.R")

#' function to simulate data under the power prior (PP) assumption
#' generate current and historical data sets based on `theta` and bootstrap samples from `X`
#'
#' @param prop.cens     desired censoring proportion
#' @param nevents       desired number of events in simulated current data
#' @param n0events      desired number of events in simulated historical data
#' @param X             design matrix w/ column names being "(intercept)", "treatment", and covariate
#'                      names (if any). Bootstrap samples from the covariates of `X` will be used as
#'                      the covariates for simulated data. Defaults to design matrix from E1690 data.
#' @param theta         a vector of scale parameter and regression coefficients, e.g., the MLEs
#'                      from fitting a Weibull AFT model on current/historical data.
#' @param trt.prob      probability of being assigned to the treated group. The treatment indicator
#'                      for simulated data will be generated using Bernoulli distribution with probability
#'                      being `trt.prob`. Defaults to 0.5.
#'
sim.PP <- function(
    prop.cens,
    nevents,
    n0events,
    X             = design.mtx.curr,
    theta         = mle.curr,
    trt.prob      = 0.5
){
  n         <- ceiling( nevents / (1 - prop.cens) )
  # simulate stacked data
  df.all    <- get.weibull.surv(
    X = X,
    theta = theta,
    prop.cens = prop.cens,
    nevents = nevents + n0events,
    trt.prob = trt.prob,
    bootstrap = TRUE
  )
  
  data      <- df.all[(1:n), ]
  histdata  <- df.all[-(1:n), ]
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}
