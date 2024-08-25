library(collapse)

#' Helper functions for computing log-likelihood matrix based on posterior samples
#'
#' Compute mean from linear predictor in a GLM
#' @param eta  linear predictor
#' @param link integer giving link function
#' @noRd
get_lp2mean = function(eta, link) {
  if (link == 1){
    return(eta) # identity link
  }else if(link == 2){
    return(exp(eta)) # log link
  }else if(link == 3){
    return(binomial('logit')$linkinv(eta)) # logit link
  }else if(link == 4){
    return(1/eta) # inverse link
  }else if(link == 5){
    return(binomial('probit')$linkinv(eta)) # probit link
  }else if(link == 6){
    return(binomial('cauchit')$linkinv(eta)) # cauchit link
  }else if(link == 7){
    return(binomial('cloglog')$linkinv(eta)) # complementary log-log link
  }else if(link == 8){
    return(eta^2) # sqrt link
  }else if(link == 9){
    return(1/sqrt(eta)) # 1/mu^2 link
  }else{
    stop("Link not supported")
  }
}

#' compute log-likelihood matrix for normal GLM
#' @param y           response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     an M x 1 vector of dispersion parameters (variance)
#' @noRd
normal_glm_loglik_mtx = function(y, beta_mtx, X, link, offs, phi_vec) {
  n         = length(y)
  theta_mtx = tcrossprod(beta_mtx, X) %r+% offs

  if ( link != 1 ){
    theta_mtx = get_lp2mean(theta_mtx, link)
  }
  sd_vec    = sqrt(phi_vec)
  res       = sapply(1:nrow(beta_mtx), function(m){
    dnorm(y, mean = theta_mtx[m, ], sd = sd_vec[m], log = T)
  })
  return( t(res) )
}

#' compute log-likelihood matrix for bernoulli GLM
#' @param y           response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     dispersion parameter (phi_vec = 1)
#' @noRd
bernoulli_glm_loglik_mtx = function(y, beta_mtx, X, link, offs, phi_vec = 1) {
  n         = length(y)
  theta_mtx = tcrossprod(beta_mtx, X) %r+% offs

  if ( link != 3 ){
    theta_mtx = binomial('logit')$linkfun( get_lp2mean(theta_mtx, link) )
  }
  #res      = sweep(theta_mtx, 2, y, "*") - log1p(exp(theta_mtx))
  res       = theta_mtx %r*% y - log1p(exp(theta_mtx))
  return(res)
}

#' compute log-likelihood matrix for poisson GLM
#' @param y           response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     dispersion parameter (phi_vec = 1)
#' @noRd
poisson_glm_loglik_mtx = function(y, beta_mtx, X, link, offs, phi_vec = 1) {
  n         = length(y)
  theta_mtx = tcrossprod(beta_mtx, X) %r+% offs

  if ( link != 2 ){
    theta_mtx = log( get_lp2mean(theta_mtx, link) )
  }
  res       = ( theta_mtx %r*% y - exp(theta_mtx) ) %r-%  lgamma(y + 1)
  return(res)
}

#' compute log-likelihood matrix for gamma GLM
#' @param y           response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     an M x 1 vector of dispersion parameters
#' @noRd
gamma_glm_loglik_mtx = function(y, beta_mtx, X, link, offs, phi_vec) {
  n         = length(y)
  theta_mtx = tcrossprod(beta_mtx, X) %r+% offs
  tau_vec   = 1 / phi_vec # shape parameter

  if ( link != 4 ){
    theta_mtx = 1 / get_lp2mean(theta_mtx, link)
  }
  rate_mtx  = theta_mtx * tau_vec
  res       = sapply(1:nrow(beta_mtx), function(m){
    dgamma(y, shape = tau_vec[m], rate = rate_mtx[m, ], log = T)
  })
  return( t(res) )
}

#' compute log-likelihood matrix for inverse-Gaussian GLM
#' @param y    response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     an M x 1 vector of dispersion parameters
#' @noRd
invgauss_glm_loglik_mtx = function(y, beta_mtx, X, link, offs, phi_vec) {
  n              = length(y)
  theta_mtx      = tcrossprod(beta_mtx, X) %r+% offs
  sqrt_theta_mtx = sqrt(theta_mtx)
  tau_vec        = 1 / phi_vec # shape parameter
  log_tau_vec    = log(tau_vec)
  log_2pi        = log(2*pi)
  log_y          = log(y)
  inv_sqrt_y     = 1/sqrt(y)

  if ( link != 9 ){
    theta_mtx = 1 / ( get_lp2mean(theta_mtx, link)^2 )
  }
  res       = sapply(1:nrow(beta_mtx), function(m){
    0.5 * (
      log_tau_vec[m] - log_2pi - 3*log_y - tau_vec[m] * ( (y * sqrt_theta_mtx[m, ] - 1) * inv_sqrt_y )^2
    )
  })
  return( t(res) )
}

#' wrapper function to compute log-likelihood matrix for a given link function and a given distribution
#' @param y           response vector
#' @param beta_mtx    an M x p matrix of regression coefficients, where M is the number of posterior samples, and p
#'                    is the number of regression coefficients per sampling iteration.
#' @param X           design matrix
#' @param dist        integer giving distribution
#' @param link        integer giving link function
#' @param offs        offset
#' @param phi_vec     an M x 1 vector of dispersion parameters
#'
glm_loglik_mtx = function(y, beta_mtx, X, dist, link, offs, phi_vec) {
  # Compute likelihood
  if (dist == 1) {     # Bernoulli
    return( bernoulli_glm_loglik_mtx(y, beta_mtx, X, link, offs, phi_vec = 1) )
  }else if (dist == 2) {  # Poisson
    return( poisson_glm_loglik_mtx(y, beta_mtx, X, link, offs, phi_vec = 1) )
  }else if (dist == 3) {  # Normal
    return( normal_glm_loglik_mtx(y, beta_mtx, X, link, offs, phi_vec) )
  }else if (dist == 4) { # Gamma
    return( gamma_glm_loglik_mtx(y, beta_mtx, X, link, offs, phi_vec) )
  }else if (dist == 5) { # Inverse-Gaussian
    return( invgauss_glm_loglik_mtx(y, beta_mtx, X, link, offs, phi_vec) )
  }else{
    stop("Distribution not supported")
  }
}
