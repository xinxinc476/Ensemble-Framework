//
//  Stratified power prior for GLMs
//  Power parameter fixed
//
functions{
  //' Compute mean from linear predictor in a GLM
  //'
  //' @param eta linear predictor
  //' @param link integer giving link function
  vector lp2mean(vector eta, int link) {
    if (link == 1) return(eta);                        // identity link
    else if (link == 2) return exp(eta);               // log link
    else if (link == 3) return inv_logit(eta);         // logit link
    else if (link == 4) return inv(eta);               // inverse link
    else if (link == 5) return Phi_approx(eta);        // probit link
    else if (link == 6) return atan(eta) / pi() + 0.5; // cauchit link
    else if (link == 7) return inv_cloglog(eta);       // complementary log-log link
    else if (link == 8) return square(eta);            // sqrt link
    else if (link == 9) return inv_sqrt(eta);          // 1/mu^2 link
    else reject("Link not supported");
    return eta; // never reached
  }

  real normal_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 1 )
      theta = lp2mean(theta, link);
    return normal_lpdf(y | theta, sqrt(phi) );
  }

  real bernoulli_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 3 )
      theta = logit( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( log1p_exp(theta) );
  }

  real poisson_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 2 )
      theta = log( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( exp(theta) + lgamma(y + 1) );
  }

  real gamma_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    real tau        = inv(phi); // shape parameter
    vector[n] theta = X * beta + offs;
    if ( link != 4 )
      theta = inv( lp2mean(theta, link) );
    return gamma_lpdf(y | tau, tau * theta );
  }

  real invgauss_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n                 = rows(y);
    real tau              = inv(phi); // shape parameter
    real log_2pi          = 1.837877066409345483560659;  // log(2*pi)
    vector[n] theta       = X * beta + offs;
    if ( link != 9 )
      theta = inv_square( lp2mean(theta, link) );
    return 0.5 * (
              n * (log(tau) - log_2pi) - 3 * sum(log(y))
            - tau * dot_self( (y .* sqrt(theta) - 1) .* inv_sqrt(y) )
          );
  }

  real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offs) {
    // Compute likelihood
    if (dist == 1) {     // Bernoulli
      return bernoulli_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 2) {  // Poisson
      return poisson_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 3) {  // Normal
      return normal_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 4) { // Gamma
      return gamma_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 5) { // Inverse-Gaussian
      return invgauss_glm_lp(y, beta, phi, X, link, offs);
    }
    else reject("Distribution not supported");
    return 0; // never reached;
  }
}

data {
  int<lower=0>                                  n;  // current data sample size
  int<lower=0>                                  p;  // number of covariates (including intercept and treatment indicator)
  vector[n]                                     y;  // current data response
  matrix[n,p]                                   X;  // current data design matrix
  vector[n]                                  offs;  // current data offset
  int<lower=0>                                 n0;  // historical data sample size
  vector[n0]                                   y0;  // historical data response
  matrix[n0,p]                                 X0;  // historical data design matrix
  vector[n0]                                offs0;  // historical data offset

  int<lower = 0>                                K;  // number of strata
  array[K] int<lower=0, upper=n>   start_idx_curr;  // starting index of each stratum in current data
  array[K] int<lower=0, upper=n>     end_idx_curr;  // ending index of each stratum in current data
  array[K] int<lower=0, upper=n0>  start_idx_hist;  // starting index of each stratum in historical data
  array[K] int<lower=0, upper=n0>    end_idx_hist;  // ending index of each stratum in historical data

  vector<lower=0, upper=1>[K]                 a0s;  // power prior parameter for each stratum

  matrix[p,K]                           beta_mean;  // mean parameter for normal prior on beta
  matrix<lower=0>[p,K]                    beta_sd;  // sd parameter for normal prior on beta
  vector[K]                             disp_mean;  // location parameter for half-normal prior on dispersion
  vector<lower=0>[K]                      disp_sd;  // scale parameter for half-normal prior on dispersion

  int<lower=1, upper=5>                      dist;
  int<lower=1, upper=9>                      link;
}

transformed data {
  real lognc_disp  = normal_lccdf(0 | disp_mean, disp_sd); // log normalizing constant for half-normal prior on dispersion
}

parameters {
  matrix[p,K]  beta; // each column is a stratum-specific vector of coefficients
  vector<lower=0>[(dist > 2) ? K :  0] dispersion;
}

model {
  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      // current data likelihood
      target += glm_lp(y[ start_idx_curr[k]:end_idx_curr[k] ],
        beta[, k], 1.0, X[ start_idx_curr[k]:end_idx_curr[k], ], dist, link,
        offs[ start_idx_curr[k]:end_idx_curr[k] ]);

      // power prior
      target += a0s[k] * glm_lp(y0[ start_idx_hist[k]:end_idx_hist[k] ],
        beta[, k], 1.0, X0[ start_idx_hist[k]:end_idx_hist[k], ], dist, link,
        offs0[ start_idx_hist[k]:end_idx_hist[k] ]);
      // initial prior for beta
      target += normal_lpdf(beta[, k] | beta_mean[, k], beta_sd[, k]);
    }
  }
  else {
    // half-normal prior for dispersion parameters
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;

    for ( k in 1:K ) {
      // current data likelihood
      target += glm_lp(y[ start_idx_curr[k]:end_idx_curr[k] ],
        beta[, k], dispersion[k], X[ start_idx_curr[k]:end_idx_curr[k], ], dist, link,
        offs[ start_idx_curr[k]:end_idx_curr[k] ]);

      // power prior
      target += a0s[k] * glm_lp(y0[ start_idx_hist[k]:end_idx_hist[k] ],
        beta[, k], dispersion[k], X0[ start_idx_hist[k]:end_idx_hist[k], ], dist, link,
        offs0[ start_idx_hist[k]:end_idx_hist[k] ]);
      // initial prior for beta
      target += normal_lpdf(beta[, k] | beta_mean[, k], beta_sd[, k]);
    }
  }
}

