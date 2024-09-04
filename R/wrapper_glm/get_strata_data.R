library(psrwe)

#' Estimate propensity scores (PS), assign subjects into different strata based on PS, and obtain stratum-specific power
#' prior parameter values (a0s) using psrwe package
#'
#' @param ps_fml_covs     right side of the propensity score (PS) formula involving the covariates, e.g., ~ sex + age.
#' @param v_arm           column name corresponding to arm (treatment v.s. control) assignment.
#' @param borrow_ctl_only whether to borrow historical information on control arm only. If FALSE, then historical information on
#'                        both the treatment and control arm are borrowed, and the strata are formed based on the PS of all
#'                        subjects in current study. Otherwise, only information on control arm will be borrowed, and the strata
#'                        will be formed based on the PS of subjects in the current study control arm. Defaults to TRUE.
#'
get.strata.data = function(
    data.list,
    ps_fml_covs,
    ps_method       = c("logistic", "randomforest"),
    v_arm           = "trt",
    ctl_arm_level   = 0,
    borrow_ctl_only = TRUE,
    nstrata         = 5,
    trim_ab         = c("both", "above", "below", "none"),
    total_borrow    = NULL,
    method          = c("distance", "inverse_distance", "n_current", "n_external")
) {
  data             = data.list[[1]]
  ## stack all historical data sets into one historical data set
  histdata         = do.call(rbind, data.list[-1])
  if ( borrow_ctl_only ){
    histdata = histdata[ histdata[, v_arm] == ctl_arm_level, ]
  }
  n              = nrow(data) ## current data sample size
  n0             = nrow(histdata) ## historical data sample size
  data.all       = as.data.frame( rbind(data, histdata) )
  data.all$study = rep(c(1, 0), times = c(n, n0))

  ## estimate propensity scores (PS) via psrwe::psrwe_est() function
  data_ps        = psrwe::psrwe_est(
    data = data.all,
    ps_fml = as.formula( paste0("study", ps_fml_covs) ),
    ps_method = ps_method,
    v_grp = "study",
    cur_grp_level = 1,
    v_arm = v_arm,
    ctl_arm_level = ctl_arm_level,
    stra_ctl_only = borrow_ctl_only,
    nstrata = nstrata,
    trim_ab = trim_ab
  )

  ## split the number of subjects to be borrowed among strata based on PS via psrwe::psrwe_borrow() function
  if ( is.null(total_borrow) )
    total_borrow = n0
  ps_bor     = psrwe::psrwe_borrow(
    dtaps = data_ps,
    total_borrow = total_borrow,
    method = method
  )
  a0.strata = ps_bor$Borrow$Alpha ## a0 for each stratum
  a0.strata[is.na(a0.strata)] = 0

  ## obtain data sets with strata assignments
  data.all.strata = ps_bor$data
  ## remove subjects in historical study who are trimmed based on PS
  data.all.strata = data.all.strata[!is.na(data.all.strata$`_strata_`), ]
  strata          = as.integer( data.all.strata$`_strata_` )

  data.list.new   = list(
    curr = data.all.strata[data.all.strata$study == 1, colnames(data)],
    hist = data.all.strata[data.all.strata$study == 0, colnames(histdata)]
  )
  strata.list     = list(
    curr = strata[ data.all.strata$study == 1 ],
    hist = strata[ data.all.strata$study == 0 ]
  )
  res             = list(
    "data.list"        = data.list.new,
    "strata.list"      = strata.list,
    "a0.strata"        = a0.strata,
    "res.psrwe.est"    = data_ps,
    "res.psrwe.borrow" = ps_bor
  )
  return(res)
}
