#' function to obtain (posterior) samples from the ensemble method based on (posterior) samples
#' from individual priors and (posterior) weights
#'
sample.ensemble = function(wts, samples.mtx){
  wts <- as.numeric(wts)
  ## draw n i.i.d. samples (c0) from categorical distribution with probability being `wts`
  c0 <- sample(x = seq_len(ncol(samples.mtx)), size = nrow(samples.mtx), replace = T,
               prob = wts)
  
  models <- unique(c0)
  res.samples <- lapply(models, function(j){
    nsample = sum(c0 == j)
    return( sample(samples.mtx[, j], size = nsample, replace = T) )
  })
  res.samples <- unlist(res.samples)
  return(res.samples)
}