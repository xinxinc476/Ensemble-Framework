library(hdbayes)

#' check if the input data are in appropriate forms for stratified power prior (PP)
#'
#' @noRd
data.checks.stratified.pp = function(
    formula, family, data.list, strata.list, offset.list
) {
  hdbayes:::data.checks(formula, family, data.list, offset.list)

  ## check strata.list
  if ( !( is.list(strata.list) ) )
    stop("strata.list must be a list of vectors")
  if ( length(strata.list) != length(data.list) )
    stop("strata.list and data.list must have equal lengths")
  for( i in seq_len( length(strata.list) ) ){
    if ( !( is.vector(strata.list[[i]]) ) )
      stop("element ", i, " in strata.list must be a vector")
    if ( any( is.na(strata.list[[i]]) ) )
      stop("element ", i, " in strata.list cannot contain missing values")
    strata.list[[i]] = as.integer( strata.list[[i]] )
    if ( any( strata.list[[i]] < 1 ) )
      stop("element ", i, " in strata.list must be a vector of integers ranging from 1 to the number of strata")
    if ( length(strata.list[[i]]) != nrow(data.list[[i]]) )
      stop("the length of element ", i, " in strata.list must be equal to the number of rows in element ", i, " in data.list")
  }
}


#' get Stan data for stratified power prior (PP)
#'
#' @noRd
get.stan.data.stratified.pp = function(
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
    disp.sd           = NULL
) {
  data.checks.stratified.pp(
    formula, family, data.list, strata.list, offset.list
  )

  ## get model information
  data             = data.list[[1]]
  data$stratum     = as.integer( strata.list[[1]] )
  offset           = offset.list[[1]]
  ## stack all historical data sets into one historical data set
  histdata         = do.call(rbind, data.list[-1])
  histdata$stratum = as.integer( unlist(strata.list[-1]) )
  hist.offset.list = offset.list[-1]

  n                = nrow(data)
  n0               = nrow(histdata)
  ## default offset for each data set is a vector of 0s
  if ( is.null(offset) ){
    offset = rep(0, n)
  }
  if ( is.null(hist.offset.list) ){
    offset0 = rep(0, n0)
  }else {
    offset0 = as.numeric( unlist(hist.offset.list) )
  }

  if( !is.sorted ){
    ## re-order data and offset by stratum
    data           = data[order(data$stratum), ]
    offset         = offset[order(data$stratum)]
    ## re-order histdata and offset0 by stratum
    histdata       = histdata[order(histdata$stratum), ]
    offset0        = offset0[order(histdata$stratum)]
  }

  y                = data[, all.vars(formula)[1]]
  X                = model.matrix(formula, data)
  p                = ncol(X)
  y0               = histdata[, all.vars(formula)[1]]
  X0               = model.matrix(formula, histdata)
  fam.indx         = hdbayes:::get.dist.link(family)
  dist             = fam.indx[1]
  link             = fam.indx[2]

  ## get the number of strata
  K = as.integer( max( data$stratum, histdata$stratum ) )

  ## get starting and ending indices for each strata
  num.obs.curr     = as.numeric( table(data$stratum) )
  end.index.curr   = cumsum(num.obs.curr)
  start.index.curr = c(1, end.index.curr[-K] + 1)
  num.obs.hist     = as.numeric( table(histdata$stratum) )
  end.index.hist   = cumsum(num.obs.hist)
  start.index.hist = c(1, end.index.hist[-K] + 1)

  ## check a0.strata values
  if ( !( is.vector(a0.strata) & (length(a0.strata) %in% c(1, K)) ) )
    stop("a0.strata must be a scalar or a vector of length ", K)
  a0.strata = hdbayes:::to.vector(param = a0.strata, len = K)
  if ( any(a0.strata < 0 | a0.strata > 1 ) )
    stop("Each element of a0.strata must be a scalar between 0 and 1")

  ## default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) ){
    beta.mean = matrix(rep(0, p*K), ncol=K)
  }else {
    if ( length(beta.mean) == 1 ){
      beta.mean = matrix(rep(as.numeric(beta.mean), p*K), ncol=K)
    }else {
      if ( !( is.matrix(beta.mean) ) )
        stop("beta.mean must be a matrix")
      if ( nrow(beta.mean) != p )
        stop("beta.mean must have ", p, " row(s) if it is not NULL")
      if ( ncol(beta.mean) != K )
        stop("beta.mean must have ", K, " column(s) if it is not NULL")
      beta.mean = matrix(as.numeric(beta.mean), nrow = p, ncol = K)
    }
  }

  if ( is.null(beta.sd) ){
    beta.sd = matrix(rep(10, p*K), ncol=K)
  }else {
    if ( length(beta.sd) == 1 ){
      beta.sd = matrix(rep(as.numeric(beta.sd), p*K), ncol=K)
    }else {
      if ( !( is.matrix(beta.sd) ) )
        stop("beta.sd must be a matrix")
      if ( nrow(beta.sd) != p )
        stop("beta.sd must have ", p, " row(s) if it is not NULL")
      if ( ncol(beta.sd) != K )
        stop("beta.sd must have ", K, " column(s) if it is not NULL")
      beta.sd = matrix(as.numeric(beta.sd), nrow = p, ncol = K)
    }
  }

  ## default half-normal hyperprior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) %in% c(1, K)) ) )
      stop("disp.mean must be a scalar or a vector of length ", K, " if disp.mean is not NULL")
  }
  disp.mean = hdbayes:::to.vector(param = disp.mean, default.value = 0, len = K)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) %in% c(1, K)) ) )
      stop("disp.sd must be a scalar or a vector of length ", K, " if disp.sd is not NULL")
  }
  disp.sd = hdbayes:::to.vector(param = disp.sd, default.value = 10, len = K)

  standat = list(
    'n'              = n,
    'p'              = p,
    'y'              = y,
    'X'              = X,
    'offs'           = offset,
    'n0'             = n0,
    'y0'             = y0,
    'X0'             = X0,
    'offs0'          = offset0,
    'K'              = K,
    'start_idx_curr' = start.index.curr,
    'end_idx_curr'   = end.index.curr,
    'start_idx_hist' = start.index.hist,
    'end_idx_hist'   = end.index.hist,
    'a0s'            = a0.strata,
    'beta_mean'      = beta.mean,
    'beta_sd'        = beta.sd,
    'disp_mean'      = disp.mean,
    'disp_sd'        = disp.sd,
    'dist'           = dist,
    'link'           = link
  )
  return(standat)
}
