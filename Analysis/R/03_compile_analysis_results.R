## Compile analysis results
library(dplyr)
library(posterior)
library(loo)

## source wrapper functions for computing the difference in 2-year RFS probability for treated v.s. untreated,
## i.e., P(T > 2 | A = 1) - P(T > 2 | A = 0)
wrapper.dir <- 'wrapper'
source(file.path(wrapper.dir, "pwe_get_trt_effect.R"))
source(file.path(wrapper.dir, "curepwe_get_trt_effect.R"))
source(file.path(wrapper.dir, "sample_ensemble.R"))

grid        <- readRDS("Analysis/R/grid.rds") %>%
  filter( arm == 0 ) %>%
  dplyr::select(!arm)
results.dir <- "Analysis/Results"
file.list   <- list.files(results.dir, pattern = ".rds")

surv.trt.list  <- list()
surv.ctl.list  <- list()
lpd.trt.list   <- list()
lpd.ctl.list   <- list()
logml.trt.vals <- vector(length = nrow(grid))
logml.ctl.vals <- vector(length = nrow(grid))
surv.diff.list <- list()

for (i in 1:nrow(grid)) {
  grid.id <- grid[i, ]
  prior   <- grid.id$prior
  model   <- grid.id$model
  
  # for IFN arm
  pattern <- paste0('_', model, "_", prior, "_arm_", 1, '_nintervals_', grid.id$J)
  res.i   <- readRDS(file.path(results.dir, file.list[ grep(pattern, file.list) ]))
  d       <- res.i$draws
  lpd.trt.list[[i]] <- res.i$res.loo$pointwise[,"elpd_loo"]
  logml.trt.vals[i] <- res.i$res.logml$logml
  
  if ( model == "PWE"){
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.pwe.psipp(
        t = 2,
        post.samples = d
      )
      
    }else{
      d.surv <- get.surv.prob.pwe(
        t = 2,
        post.samples = d
      )
    }
    
  }else{
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.curepwe.psipp(
        t = 2,
        post.samples = d
      )
      
    }else{
      d.surv <- get.surv.prob.curepwe(
        t = 2,
        post.samples = d
      )
      
    }
  }
  surv.trt.list[[i]] <- d.surv
  
  # for Control arm
  pattern <- paste0('_', model, "_", prior, "_arm_", 0, '_nintervals_', grid.id$J)
  res.i   <- readRDS(file.path(results.dir, file.list[ grep(pattern, file.list) ]))
  d       <- res.i$draws
  lpd.ctl.list[[i]] <- res.i$res.loo$pointwise[,"elpd_loo"]
  logml.ctl.vals[i] <- res.i$res.logml$logml
  
  if ( model == "PWE"){
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.pwe.psipp(
        t = 2,
        post.samples = d
      )
      
    }else{
      d.surv <- get.surv.prob.pwe(
        t = 2,
        post.samples = d
      )
    }
    
  }else{
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.curepwe.psipp(
        t = 2,
        post.samples = d
      )
      
    }else{
      d.surv <- get.surv.prob.curepwe(
        t = 2,
        post.samples = d
      )
      
    }
  }
  surv.ctl.list[[i]] <- d.surv
  
  d.surv.diff         <- surv.trt.list[[i]] - surv.ctl.list[[i]]
  surv.diff.list[[i]] <- d.surv.diff
  estim.i <- round(c(mean = mean(d.surv.diff), sd = sd(d.surv.diff),
                     quantile2(d.surv.diff, probs = c(0.5, 0.025, 0.975)),
                     prob_greater_0 = mean(d.surv.diff > 0)), 3)
  prior.i <- prior
  model.i <- paste0(model, " (J = ", grid.id$J, ")")
  
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
}
res <- list(
  surv.trt.list  = surv.trt.list,
  surv.ctl.list  = surv.ctl.list,
  lpd.trt.list   = lpd.trt.list,
  lpd.ctl.list   = lpd.ctl.list,
  logml.trt.vals = logml.trt.vals,
  logml.ctl.vals = logml.ctl.vals,
  surv.diff.list = surv.diff.list,
  estim = estim,
  priors = priors,
  models = models
)

## compute the 2-year RFS probabilities for IFN arm using the ensemble methods
## assume uniform initial construct probabilities

## weights from BMA, Pseudo-BMA, Pseudo-BMA+, and Stacking
post.prob    <- as.numeric( bridgesampling::post_prob(res$logml.trt.vals) )
lpd_point    <- do.call(cbind, res$lpd.trt.list)
pbma_wts     <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts  <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)

surv.trt.mtx      <- do.call(cbind, res$surv.trt.list)
surv.trt.bma      <- sample.ensemble(wts = post.prob, samples.mtx = surv.trt.mtx)
surv.trt.pbma     <- sample.ensemble(wts = pbma_wts, samples.mtx = surv.trt.mtx)
surv.trt.pbma.BB  <- sample.ensemble(wts = pbma_BB_wts, samples.mtx = surv.trt.mtx)
surv.trt.stacking <- sample.ensemble(wts = stacking_wts, samples.mtx = surv.trt.mtx)
surv.trt.list.all <- append(res$surv.trt.list,
                            list(bma = surv.trt.bma, pbma = surv.trt.pbma,
                                 pbmaBB = surv.trt.pbma.BB, stacking = surv.trt.stacking))

wts.tab.trt  <- data.frame(
  model = res$models,
  prior = res$priors,
  "BMA" = post.prob,
  "Pseudo-BMA" = pbma_wts,
  "Pseudo-BMA+" = pbma_BB_wts, ## Pseudo-BMA with Bayesian bootstrap
  "Stacking" = stacking_wts
)
rownames(wts.tab.trt) <- NULL
colnames(wts.tab.trt)[4:5] <- c("Pseudo-BMA", "Pseudo-BMA+")

model_order <- c(paste0("CurePWE (J = ", 2:9, ")"), paste0("PWE (J = ", 2:9, ")"))

wts.tab.trt <- wts.tab.trt %>%
  mutate(model = factor(model, levels = model_order)) %>%
  arrange(prior, model)


## compute the 2-year RFS probabilities for Control arm using the ensemble methods
## assume uniform initial construct probabilities

## weights from BMA, Pseudo-BMA, Pseudo-BMA+, and Stacking
post.prob    <- as.numeric( bridgesampling::post_prob(res$logml.ctl.vals) )
lpd_point    <- do.call(cbind, res$lpd.ctl.list)
pbma_wts     <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts  <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)

surv.ctl.mtx      <- do.call(cbind, res$surv.ctl.list)
surv.ctl.bma      <- sample.ensemble(wts = post.prob, samples.mtx = surv.ctl.mtx)
surv.ctl.pbma     <- sample.ensemble(wts = pbma_wts, samples.mtx = surv.ctl.mtx)
surv.ctl.pbma.BB  <- sample.ensemble(wts = pbma_BB_wts, samples.mtx = surv.ctl.mtx)
surv.ctl.stacking <- sample.ensemble(wts = stacking_wts, samples.mtx = surv.ctl.mtx)
surv.ctl.list.all <- append(res$surv.ctl.list,
                            list(bma = surv.ctl.bma, pbma = surv.ctl.pbma,
                                 pbmaBB = surv.ctl.pbma.BB, stacking = surv.ctl.stacking))

wts.tab.ctl  <- data.frame(
  model = res$models,
  prior = res$priors,
  "BMA" = post.prob,
  "Pseudo-BMA" = pbma_wts,
  "Pseudo-BMA+" = pbma_BB_wts, ## Pseudo-BMA with Bayesian bootstrap
  "Stacking" = stacking_wts
)
rownames(wts.tab.ctl) <- NULL
colnames(wts.tab.ctl)[4:5] <- c("Pseudo-BMA", "Pseudo-BMA+")

model_order <- c(paste0("CurePWE (J = ", 2:9, ")"), paste0("PWE (J = ", 2:9, ")"))

wts.tab.ctl <- wts.tab.ctl %>%
  mutate(model = factor(model, levels = model_order)) %>%
  arrange(prior, model)

## compute the difference in 2-year RFS probabilities using the ensemble methods
d.surv.diff.bma      <- surv.trt.bma - surv.ctl.bma
d.surv.diff.pbma     <- surv.trt.pbma - surv.ctl.pbma
d.surv.diff.pbma.BB  <- surv.trt.pbma.BB - surv.ctl.pbma.BB
d.surv.diff.stacking <- surv.trt.stacking - surv.ctl.stacking
surv.diff.list.all   <- append(res$surv.diff.list,
                               list(bma = d.surv.diff.bma, pbma = d.surv.diff.pbma,
                                    pbmaBB = d.surv.diff.pbma.BB, stacking = d.surv.diff.stacking))

estim2 <- lapply(list(d.surv.diff.bma, d.surv.diff.pbma, d.surv.diff.pbma.BB, d.surv.diff.stacking), function(l){
  round(c(mean = mean(l), sd = sd(l), quantile2(l, probs = c(0.5, 0.025, 0.975)),
          prob_greater_0 = mean(l > 0)), 3)
})
estim2 <- do.call(rbind, estim2) %>%
  as.data.frame()

estim.all <- rbind(res$estim, estim2) %>%
  as.data.frame()
estim.all <- cbind(model = c(res$models, rep("Ensemble", 4)),
                   prior = c(res$priors, c("BMA", "Pseudo-BMA", "Pseudo-BMA+", "Stacking")),
                   estim.all) %>%
  as.data.frame()
rownames(estim.all) <- NULL


## save compiled analysis results
res.all <- list(
  surv.trt.list.all  = surv.trt.list.all,
  surv.ctl.list.all  = surv.ctl.list.all,
  surv.diff.list.all = surv.diff.list.all,
  wts.tab.trt = wts.tab.trt,
  wts.tab.ctl = wts.tab.ctl,
  estim.all = estim.all,
  models = res$models,
  priors = res$priors
)
saveRDS(res.all, "Analysis/Results/compiled_analysis_results.rds")
