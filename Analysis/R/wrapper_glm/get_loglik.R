library(hdbayes)
source("Analysis/R/wrapper_glm/glm_loglik_mtx.R")

#' Pointwise log-likelihood matrix of a GLM under stratified power prior (PP)
#'
#' @description This function computes the pointwise log-likelihood matrix of a GLM based on posterior samples under the
#' stratified power prior (PP). The output can be used to efficiently approximate leave-one-out cross-validation (LOO)
#' via [loo::loo()].
#'
#' @param post.samples      output from [glm.stratified.pp()] giving posterior samples of a GLM under the stratified power
#'                          prior (PP), with an attribute called 'data' which includes the list of variables specified in
#'                          the data block of the Stan program.
#'
#' @return
#'  The function returns a \eqn{M \times n} log-likelihood matrix, where \eqn{M} is the number of posterior samples, and
#'  \eqn{n} is the sample size of current data.
#'
get.loglik.glm.stratified.pp = function(post.samples){
  stan.data = attr(post.samples, 'data')
  n         = stan.data$n
  p         = stan.data$p
  K         = stan.data$K
  y         = stan.data$y
  X         = stan.data$X
  offs      = stan.data$offs
  start.idx = stan.data$start_idx_curr
  end.idx   = stan.data$end_idx_curr
  dist      = stan.data$dist
  link      = stan.data$link

  beta.names = paste0( colnames(X), '_stratum_', rep(1:K, each = p) )
  beta       = suppressWarnings(
    as.matrix( post.samples[, beta.names, drop=F] )
  )

  if ( dist > 2 ){
    disp.names = paste0( 'dispersion', '_stratum_', 1:K )
    dispersion = suppressWarnings(
      as.matrix( post.samples[, disp.names, drop=F] )
    )

    loglik = lapply(1:K, function(k){
      y.stratum    = y[ start.idx[k]:end.idx[k] ]
      X.stratum    = X[ start.idx[k]:end.idx[k], ]
      offs.stratum = offs[ start.idx[k]:end.idx[k] ]
      disp.stratum = dispersion[, k]
      beta.stratum = beta[, paste0( colnames(X), '_stratum_', k)]
      glm_loglik_mtx(
        y = y.stratum, beta_mtx = beta.stratum, X = X.stratum,
        dist = dist, link = link, offs = offs.stratum, phi_vec = disp.stratum
      )
    })
  }else {
    loglik = lapply(1:K, function(k){
      y.stratum    = y[ start.idx[k]:end.idx[k] ]
      X.stratum    = X[ start.idx[k]:end.idx[k], ]
      offs.stratum = offs[ start.idx[k]:end.idx[k] ]
      beta.stratum = beta[, paste0( colnames(X), '_stratum_', k)]
      glm_loglik_mtx(
        y = y.stratum, beta_mtx = beta.stratum, X = X.stratum,
        dist = dist, link = link, offs = offs.stratum, phi_vec = 1
      )
    })
  }
  loglik = do.call(cbind, loglik)
  return(loglik)
}


#' Pointwise log-likelihood matrix of a GLM under various priors, including power prior (PP), normal/half-normal prior,
#' normalized power prior(NPP), normalized asymptotic power prior (NAPP), commensurate prior, Bayesian hierarchical
#' model (BHM), robust meta-analytic predictive prior (RMAP), and latent exchangeability prior (LEAP)
#'
#' @description This function computes the pointwise log-likelihood matrix of a GLM based on posterior samples under the
#' various priors. The output can be used to efficiently approximate leave-one-out cross-validation (LOO) via [loo::loo()].
#'
#' @param post.samples      output from functions like [glm.pp()] in hdbayes giving posterior samples of a GLM under
#'                          different priors, with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
#'
#' @return
#'  The function returns a \eqn{M \times n} log-likelihood matrix, where \eqn{M} is the number of posterior samples, and
#'  \eqn{n} is the sample size of current data.
#'
get.loglik.glm = function(post.samples){
  stan.data = attr(post.samples, 'data')
  if ( "n" %in% names(stan.data) ){
    n         = stan.data$n
    p         = stan.data$p
    y         = stan.data$y
    X         = stan.data$X
    offs      = stan.data$offs

    if ( "matrix" %in% class(offs) ){
      offs    = offs[, 1]
    }
  }else{
    n         = stan.data$end_idx[1] ## current data sample size
    p         = stan.data$p
    y         = stan.data$y[1:n]
    X         = stan.data$X[1:n, ]
    offs      = stan.data$offs[1:n]
  }
  dist      = stan.data$dist
  link      = stan.data$link

  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )

  if ( dist > 2 ){
    disp.name  = 'dispersion'
    if ( !disp.name %in% colnames(post.samples) ){
      disp.name = "dispersion[1]"
    }
    dispersion = suppressWarnings(
      as.matrix( post.samples[, disp.name, drop=F] )
    )
    loglik = glm_loglik_mtx(
      y = y, beta_mtx = beta, X = X, dist = dist,
      link = link, offs = offs, phi_vec = dispersion
    )
  }else {
    loglik = glm_loglik_mtx(
      y = y, beta_mtx = beta, X = X, dist = dist,
      link = link, offs = offs, phi_vec = 1
    )
  }
  return(loglik)
}
