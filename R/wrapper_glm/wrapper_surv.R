## wrapper functions for fitting a Poisson GLM for survival data
## equivalent to fitting a proportional hazards model with a piece-wise constant baseline hazards

library(survival)
library(dplyr)

source("R/wrapper_glm/get_loglik.R")

#' transform survival time data to 'counting process' format
#' @param data            a data.frame
#' @param failtime.name   variable name of failure/event time
#' @param failind.name    variable name of failure/event indicator
#' @param covariate.name  a vector of variable names (other than `failtime.name` and `failind.name`) that should be retained.
#'                        If `covariate.name` is set to NULL, then all variables in the data will be retained.
#' @param breaks          cut points of time intervals
get.counting.process.data <- function(
    data,
    breaks,
    failtime.name  = "failtime",
    failind.name   = "failcens",
    covariate.name = NULL
) {
  if( is.null(covariate.name) ){
    covariate.name = colnames(data)
    covariate.name = covariate.name[ !covariate.name %in% c(failtime.name, failind.name)]
  }
  split.formula = paste(covariate.name, collapse = "+")
  split.formula <- as.formula(
    paste("Surv(", failtime.name, ",", failind.name, ")", "~", split.formula, sep = "")
  )

  ## create pseudo-observations for each combination of subjects and intervals
  data <- survival::survSplit(formula = split.formula, data = data, cut = breaks,
                              episode = "interval", start = "start")

  breaks.labels <- paste("(", c(0, round(breaks, 2)), ", ", c(round(breaks, 2), "inf"), "]", sep="")
  ninterval     <- length(unique(data$interval))
  data          <- data %>%
    mutate(
      exposure    = .data[[failtime.name]] - start,
      log_exposue = log(exposure),
      interval    = factor(interval, labels = breaks.labels[1:ninterval])
    )
  return(data)
}


#' compute DIC for GLMs under various priors
#' @param post.samples      output from functions in hdbayes like [glm.pp()] or [glm.stratified.pp()] giving posterior
#'                          samples of a GLM under different priors, with an attribute called 'data' which includes
#'                          the list of variables specified in the data block of the Stan program.
#' @param is.stratified.pp whether `post.samples` are obtained using the stratified power prior
#'
compute.DIC = function(post.samples, is.stratified.pp){
  post.samples <- rbind(colMeans(post.samples), post.samples)

  if( is.stratified.pp ){
    loglik <- get.loglik.glm.stratified.pp(post.samples)
  }else{
    loglik <- get.loglik.glm(post.samples)
  }
  loglik.draws      <- rowSums( loglik )
  loglik.postmean   <- loglik.draws[1]
  mean.loglik.draws <- mean(loglik.draws[-1])
  d.postmean        <- -2 * loglik.postmean    ## deviance evaluated at posterior mean (deviance of average)
  d.draws           <- -2 * mean.loglik.draws  ## average deviance for each parameter (average deviance)
  pd                <- d.draws - d.postmean    ## pd = num eff params = average deviance - deviance of average
  dic               <- pd + d.draws            ## DIC = nubmer eff params + average deviance
  return(dic)
}
