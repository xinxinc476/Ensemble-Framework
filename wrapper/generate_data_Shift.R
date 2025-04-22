library(survival)
library(dplyr)

source("wrapper/wrapper_sim.R")

#' function to simulate data under the parameter shift case
#' generate current data set using `theta` and bootstrap samples from `X`
#' generate historical data set using theta0 ~ N(theta, theta.sd) and bootstrap samples from `X`
#'
#' @param prop.cens     desired censoring proportion
#' @param nevents       desired number of events in simulated current data
#' @param n0events      desired number of events in simulated historical data
#' @param X             design matrix w/ column names being "(intercept)", "treatment", and covariate
#'                      names (if any). Bootstrap samples from the covariates of `X` will be used as
#'                      the covariates for simulated data. Defaults to design matrix from E1690 data.
#' @param theta         a vector of scale parameter and regression coefficients, e.g., the MLEs
#'                      from fitting a Weibull AFT model on current/historical data.
#' @param theta.sd      standard deviation of the normal distribution from which the MLEs of log(scale)
#'                      and regression coefficients used to simulate the historical data are drawn.
#'                      Defaults to 0.2.
#' @param trt.prob      probability of being assigned to the treated group. The treatment indicator
#'                      for simulated data will be generated using Bernoulli distribution with probability
#'                      being `trt.prob`. Defaults to 0.5.
#'
sim.BHM <- function(
    prop.cens,
    nevents,
    n0events,
    X             = design.mtx.curr,
    theta         = mle.curr,
    theta.sd      = 0.2,
    trt.prob      = 0.5
){
  data      <- get.weibull.surv(
    X = X,
    theta = theta,
    prop.cens = prop.cens,
    nevents = nevents,
    trt.prob = trt.prob,
    bootstrap = TRUE
  )

  theta0    <- rnorm(n = length(theta), mean = c(log(theta[1]), theta[-1]), sd = theta.sd)
  theta0[1] <- exp(theta0[1])
  histdata  <- get.weibull.surv(
    X = X,
    theta = theta0,
    prop.cens = prop.cens,
    nevents = n0events,
    trt.prob = trt.prob,
    bootstrap = TRUE
  )
  data.list <- list(curr = data, hist = histdata)
  return(data.list)
}
