## Compile analysis results
library(dplyr)
library(posterior)

results.dir <- "Analysis/Results"
file.list   <- list.files(results.dir, pattern = ".rds")

## source wrapper functions for computing the difference in relapse-free survival probability at 2 years for
## treated v.s. untreated, i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)
wrapper.dir <- 'wrapper'
source(file.path(wrapper.dir, "pwe_get_trt_effect.R"))
source(file.path(wrapper.dir, "curepwe_get_trt_effect.R"))

surv.diff.list <- list(length = length(file.list))
loo.list       <- list(length = length(file.list))
logml.list     <- list(length = length(file.list))

for (i in 1:length(file.list)) {
  res.i         <- readRDS(file.path(results.dir, file.list[i]))
  scen.i        <- res.i$scen
  loo.list[[i]] <- res.i$res.loo
  logml.list[[i]] <- res.i$res.logml
  
  d     <- res.i$draws
  prior <- scen.i$priors
  
  if ( model == "PWE"){
    if ( prior == "psipp" ){
      d.surv.diff <- get.surv.diff.2yr.pwe.psipp(
        post.samples = d
      )
      
    }else{
      d.surv.diff <- get.surv.diff.2yr.pwe(
        post.samples = d
      )
    }
    
  }else{
    if ( prior == "psipp" ){
      d.surv.diff <- get.surv.diff.2yr.curepwe.psipp(
        post.samples = d
      )
      
    }else{
      d.surv.diff <- get.surv.diff.2yr.curepwe(
        post.samples = d
      )
      
    }
  }
  surv.diff.list[[i]] <- d.surv.diff
  
  estim.i <- round(c(mean = mean(d.surv.diff), sd = sd(d.surv.diff),
                     quantile2(d.surv.diff, probs = c(0.5, 0.025, 0.975))), 3)
  prior.i <- prior
  model.i <- paste0(scen.i$models, " (J = ", scen.i$J, ")")
  
  if( i == 1 ){
    estim  <- estim.i
    priors <- prior.i
    models <- model.i
  }else{
    estim  <- rbind(estim, estim.i)
    priors <- c(priors, prior.i)
    models <- c(models, model.i)
  }
  print( paste0("######################## Completed iteration ", i, " #########################"))
  
  res <- list(
    loo.list = loo.list,
    surv.diff.list = surv.diff.list,
    logml.list = logml.list,
    estim = estim,
    priors = priors,
    models = models
  )
  saveRDS(res, file.path(results.dir, "compiled_analysis_results.rds"))
}
