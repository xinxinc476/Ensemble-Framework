## compute and plot DIC/LOOIC/log marginal likelihood
## choose the optimal number of intervals for each prior based on log marginal likelihood (maximum)
library(tidyverse)
library(ggplot2)
library(hdbayes)

res.dir   <- 'results/PWE_choose_nintervals_RFS'
file.list <- list.files(res.dir, pattern = '.rds')

priors    <- sapply(file.list, function(x){
  prior = strsplit(x, "_")[[1]][4]
  if ( prior == "pp" ){
    prior = paste(prior, strsplit(x, "_")[[1]][5], sep = "_")
  }
  return(prior)
})

nburnin  <- 8000
nsamples <- 28000
nchains  <- 1
ncores   <- 10

wrapper.dir <- '~/projects/super_prior/Analysis/R/wrapper_glm'
source(file.path(wrapper.dir, 'glm_logml_stratified_pp.R'))

res     <- readRDS(file = "~/projects/super_prior/Analysis/R/02_compiled_results_temp.rds")
i.start <- res$i + 1
res     <- res$res

for( i in i.start:length(file.list) ){
  sim.i <- readRDS(file.path(res.dir, file.list[i]))
  if ( priors[i] == "bhm") {
    logml.i <- glm.logml.map(
      post.samples = sim.i$draws,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] == "commensurate" ){
    logml.i <- glm.logml.commensurate(
      post.samples = sim.i$draws,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] == "leap" ){
    logml.i <- glm.logml.leap(
      post.samples = sim.i$draws,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] == "npp" ){
    logml.i <- glm.logml.npp(
      post.samples = sim.i$draws,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] %in% c("pp_0.2", "pp_0.5", "pp_0.8") ){
    logml.i <- glm.logml.pp(
      post.samples = sim.i$draws,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] == "psipp" ){
    logml.i <- glm.logml.stratified.pp(
      post.samples = sim.i$draws,
      iter_warmup = nburnin,
      iter_sampling = nsamples,
      chains = nchains,
      bridge.args = list(cores = ncores)
    )
  }else if ( priors[i] == "ref" ){
    logml.i <- glm.logml.post(
      post.samples = sim.i$draws,
      bridge.args = list(cores = ncores)
    )
  }

  res.i <- cbind(sim.i$scen, DIC = sim.i$dic, LOOIC = sim.i$looic, logml = logml.i$logml)

  if ( i == 1 ){
    res <- res.i
  }else{
    res <- rbind(res, res.i)
  }

  temp = list(res = res, i = i)
  saveRDS(temp, "~/projects/super_prior/Analysis/R/02_compiled_results_temp.rds")
}
saveRDS(res, "~/projects/super_prior/Analysis/R/02_compiled_results.rds")

## recode priors
res$priors <- recode(
  res$priors, ref = "Reference", bhm = "BHM", commensurate = "Commensurate",
  npp = "NPP", leap = "LEAP", `pp_0.2` = "PP w/ a0 = 0.2", `pp_0.5` = "PP w/ a0 = 0.5",
  `pp_0.8` = "PP w/ a0 = 0.8", psipp = "PSIPP"
)

res.long <- res %>%
  pivot_longer(
    cols = c("DIC", "LOOIC"),
    names_to = "metric",
    values_to = "value"
  )

res.long %>%
  ggplot(aes(x = J, y = value, linetype = metric)) +
  geom_point() +
  geom_line() +
  labs(x = 'J', linetype = 'Information Criteria') +
  facet_wrap(~priors) +
  labs(title = "DIC/LOOIC v.s. Number of Intervals (J)")

plot_DIC_LOOIC <- function(prior){
  filenames  = file.list[priors == prior]
  dic        = vector(length = length(filenames))
  looic      = vector(length = length(filenames))
  nintervals = vector(length = length(filenames))

  for (i in 1:length(filenames)){
    res      = readRDS(file.path(res.dir, filenames[i]))
    dic[i]   = res$dic
    looic[i] = res$looic
    nintervals[i] = res$scen$J
  }
  df         = data.frame(metric = c(dic, looic),
                          type   = rep(c("DIC", "LOOIC"), each = length(dic)),
                          nintervals = rep(nintervals, 2))
 p = df %>%
   ggplot(aes(x = nintervals, y = metric, linetype = type)) +
   geom_point()  +
   geom_line() +
   labs(y = 'Info Criteria', x = 'J', linetype = '',
        title = paste0("DIC/LOOIC Values versus Number of Intervals (J) for ", prior))
 return(p)

}

plt.list = lapply(unique(priors), function(x){
  plot_DIC_LOOIC(prior = x)
})

idx <- c(1, 3, 4, 5, 6, 7, 8, 9, 10, 2)
file.list <- file.list[idx]
dic <- vector(length = length(file.list))

for (i in 1:length(file.list)){
  res <- readRDS(file.path(res.dir, file.list[i]))
  dic[i] <- res$dic
}
dic

library(ggplot2)
library(scales)
df <- data.frame(J = as.integer(2:10), dic = dic[2:10])
p1 = df %>%
  ggplot(aes(x = J, y = dic)) +
  #geom_point()  +
  geom_line() +
  labs(y = 'DIC', x = 'J',
       xlim = c(2, 10), title = "DIC Values from the Proportional Hazards Regression Model") +
  scale_x_continuous(breaks = seq(2, 10, by = 2))
# width: 620, height: 400

# curePWE:
res.dir <- '/work/users/x/i/xinxinc/strapp_paper2/curePWE_choose_J'
file.list   <- list.files(res.dir, pattern = '.rds')
data_type <- sapply(file.list, function(x){
  strsplit(x, '_')[[1]][2]
})
files_normal <- file.list[which(data_type == 'normal')]
files_logistic <- file.list[which(data_type == 'logistic')]

idx <- c(7, 8, 9, 1:6)
files_normal <- files_normal[idx]
idx <- c(6:9, 1:5)
files_logistic <- files_logistic[idx]

dic_normal = dic_logistic = vector(length = length(files_normal))

for (i in 1:length(files_normal)){
  res_normal <- readRDS(file.path(res.dir, files_normal[i]))
  dic_normal[i] <- res_normal$dic
  res_logistic <- readRDS(file.path(res.dir, files_logistic[i]))
  dic_logistic[i] <- res_logistic$dic
}
dic_normal
dic_logistic

library(ggplot2)
library(scales)
df <- data.frame(J = rep(as.integer(2:10), 2), dic = c(dic_normal, dic_logistic),
                 hist.type = rep(c('normal', 'logistic'), each = length(dic_normal)))

p2=df %>% filter(hist.type == "normal") %>%
  ggplot(aes(x = J, y = dic)) +
  #geom_point()  +
  geom_line() +
  labs(y = 'DIC', x = 'J',
       xlim = c(2, 10), title = "DIC Values from the Mixture Cure Rate Model") +
  scale_x_continuous(breaks = seq(2, 10, by = 2))
# width: 620, height: 400

plt_e1690=df %>% filter(hist.type == "logistic") %>%
  ggplot(aes(x = J, y = dic)) +
  #geom_point()  +
  geom_line() +
  labs(y = 'DIC', x = 'J',
       xlim = c(2, 10), title = "DIC Values from the Mixture Cure Rate Model") +
  scale_x_continuous(breaks = seq(2, 10, by = 2))
# width: 620, height: 400

library(ggpubr)
library(gridExtra)
plt_e1694 = ggarrange(p1, p2, nrow=1, ncol=2, labels = c('(a)', '(b)'))
annotate_figure(plt_e1694, top = text_grob("For the ECOG E1694 Trial", face = "bold", size = 14))
# width: 1200, height: 400

annotate_figure(plt_e1690, top = text_grob("For the ECOG E1690 Trial", face = "bold", size = 14))
# width: 620, height: 400
