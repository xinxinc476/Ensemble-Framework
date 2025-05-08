library(tidyverse)
library(scoringutils)

#sim.dir   <- "/work/users/x/i/xinxinc/super_prior/Sims/case_Exch"
#sim.dir   <- "/work/users/x/i/xinxinc/super_prior/Sims/case_Shift"
sim.dir   <- "/work/users/x/i/xinxinc/super_prior/Sims/case_Unexch"

file.list <- list.files(path = sim.dir, pattern = '.rds')
save.dir  <- "Sims/Results"

## specify true value of treatment effect (difference in RFS probability at 2 years for treated v.s. untreated)
param_true_list <- readRDS(file = "Sims/R/trteff_true.rds")
param_true      <- param_true_list$trteff.true
rm(param_true_list)

for(i in seq_len(length(file.list))){
  res.i <- readRDS(file.path(sim.dir, file.list[i]))
  simres.trt.i <- res.i$simres.trt
  simres.wts.i <- res.i$simres.wts
  scen.i       <- res.i$scen

  simres.trt.i$param_true <- param_true
  ## compute interval score
  ## interval score = 0.05/2 * ( (q97.5-q2.5) + 2/0.05 * max(0, q2.5 - param_true) + 2/0.05 * max(0, param_true - q97.5) )
  simres.trt.i <- simres.trt.i %>%
    mutate(
      interval_score = scoringutils::interval_score(
        true_values = param_true, lower = `q2.5`, upper = `q97.5`,
        interval_range = 95, weigh = F
      ),
      diff = mean - param_true,
      diff_sq = diff^2,
      ci.ind  = ifelse(param_true >= `q2.5` & param_true <= `q97.5`, 1, 0),
      log_var = ifelse(sd > 0, log(sd^2), NA),
      ci.width = `q97.5` - `q2.5`
    )

  ## add simulation scenario
  simres.trt.i$nevents   <- scen.i$nevents
  simres.trt.i$cens.prop <- scen.i$cens.prop
  simres.trt.i$case      <- scen.i$case

  simres.wts.i$nevents   <- scen.i$nevents
  simres.wts.i$cens.prop <- scen.i$cens.prop
  simres.wts.i$case      <- scen.i$case

  if( i == 1 ){
    simres.trt <- simres.trt.i
    simres.wts <- simres.wts.i
  }else{
    simres.trt <- rbind(simres.trt, simres.trt.i)
    simres.wts <- rbind(simres.wts, simres.wts.i)
  }
}

lst <- list(
  simres.trt = simres.trt,
  simres.wts = simres.wts
)
#saveRDS(lst, file = file.path(save.dir, "results_case_PP.rds"))
#saveRDS(lst, file = file.path(save.dir, "results_case_BHM.rds"))
saveRDS(lst, file = file.path(save.dir, "results_case_UNEXCH.rds"))


rm(list = ls())
save.dir  <- "Sims/Results"
# combine results from all three cases
res_PP     <- readRDS(file = file.path(save.dir, "results_case_PP.rds"))
res_BHM    <- readRDS(file = file.path(save.dir, "results_case_BHM.rds"))
res_UNEXCH <- readRDS(file = file.path(save.dir, "results_case_UNEXCH.rds"))

lst.combined <- list(
  simres.trt = rbind(res_PP$simres.trt, res_BHM$simres.trt, res_UNEXCH$simres.trt),
  simres.wts = rbind(res_PP$simres.wts, res_BHM$simres.wts, res_UNEXCH$simres.wts)
)

res_avg_wts <- lst.combined$simres.wts %>%
  group_by(method, nevents, `cens.prop`, case) %>%
  summarise(count = n(), mean_post_prob = mean(post.prob),
            mean_pbma_wts = mean(pbma_wts),
            mean_pbma_BB_wts = mean(pbma_BB_wts),
            mean_stacking_wts = mean(stacking_wts))

res_avg_trt <- lst.combined$simres.trt %>%
  group_by(method, nevents, `cens.prop`, case)%>%
  summarise(
    count = n(),
    avg_log_var = mean(log_var)
    ,      bias = mean(diff)
    ,   log_mse = log(mean(diff_sq))
    ,    ci.cov = mean(ci.ind)
    , avg_interval_score = mean(interval_score)
    ,       avg_ci_width = mean(ci.width)
  )
res <- list(
  res_avg_wts = res_avg_wts,
  res_avg_trt = res_avg_trt
)
saveRDS(res, file = file.path(save.dir, "compiled_Sims.rds"))
