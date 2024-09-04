library(hdbayes)
library(cmdstanr)
source("R/wrapper_glm/get_stan_data.R")

#' Sample from the posterior distribution of a GLM using the stratified power prior (PP).
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For sstratified PP implementation, all historical data sets
#'                          will be stacked into one historical data set.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of `offset.list` is equal to
#'                          the length of `data.list`. The length of each element of `offset.list` is equal to the number
#'                          of rows in the corresponding element of `data.list`. Defaults to a list of vectors of 0s.
#' @param strata.list       a list of vectors giving the strata assignments for each data. The length of `strata.list` is
#'                          equal to the length of `data.list`. Each element of `strata.list` is a vector of integers
#'                          (ranging from 1 to the number of strata) with length being equal to the number of rows in the
#'                          corresponding element of `data.list`.
#' @param is.sorted         whether each element of `data.list` (and each element of `offset.list`, if not NULL) has been
#'                          sorted in ascending order based on the corresponding element in `strata.list`. Defaults to FALSE.
#' @param beta.mean         a scalar or a `p x K` matrix of mean parameters for initial priors on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept), and `K` is the number
#'                          of strata. If a scalar is provided, `beta.mean` will be a matrix of repeated elements of the
#'                          given scalar. Defaults to a matrix of 0s.
#' @param beta.sd           a scalar or a `p x K` matrix of sd parameters for initial priors on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept), and `K` is the number
#'                          of strata. If a scalar is provided, same as for `beta.mean`. Defaults to a matrix of 10s.
#' @param a0.strata         a scalar between 0 and 1 or a vector whose dimension is equal to the number of strata giving the
#'                          (fixed) power prior parameter for each stratum. Every element of the vector should be between 0
#'                          and 1. If a scalar is provided, `a0.strata` will be a vector of repeated elements of the given
#'                          scalar.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of strata giving the location
#'                          parameters for the half-normal priors on the dispersion parameters. If a scalar is provided,
#'                          same as for `a0.strata`. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of strata giving the scale parameters
#'                          for the half-normal priors on the dispersion parameters. If a scalar is provided, same as for
#'                          `a0.strata`. Defaults to a vector of 10s.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
glm.stratified.pp = function(
    formula,
    family,
    data.list,
    strata.list,
    a0.strata,
    offset.list       = NULL,
    is.sorted         = FALSE,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  if ( length(data.list) == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for stratified power prior (PP)
  standat = get.stan.data.stratified.pp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    strata.list    = strata.list,
    a0.strata      = a0.strata,
    offset.list    = offset.list,
    is.sorted      = is.sorted,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd
  )

  glm_stratified_pp = cmdstanr::cmdstan_model("Stan/glm_stratified_pp.stan")
  ## fit model in cmdstanr
  fit = glm_stratified_pp$sample(
    data = standat,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
    ...
  )

  ## rename parameters
  p        = standat$p
  X        = standat$X
  K        = standat$K
  oldnames = paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  newnames = paste0( colnames(X), '_stratum_', rep(1:K, each = p) )

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, paste0( 'dispersion', '_stratum_', 1:K ))
  }
  d = hdbayes:::rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
