## compute and plot DIC/LOOIC/log marginal likelihood
## choose the optimal number of intervals for each prior based on log marginal likelihood (maximum)
library(tidyverse)
library(ggplot2)
library(hdbayes)
library(latex2exp)

res.dir   <- 'Analysis/Results'
file.list <- list.files(res.dir, pattern = '.rds')

priors    <- sapply(file.list, function(x){
  prior = strsplit(x, "_")[[1]][4]
  if ( prior == "pp" ){
    prior = paste(prior, strsplit(x, "_")[[1]][5], sep = "_")
  }
  return(prior)
})

for( i in 1:length(file.list) ){
  sim.i <- readRDS(file.path(res.dir, file.list[i]))
  res.i <- cbind(sim.i$scen, DIC = sim.i$dic, LOOIC = sim.i$looic, logml = sim.i$logml)

  if ( i == 1 ){
    res <- res.i
  }else{
    res <- rbind(res, res.i)
  }
}
saveRDS(res, "Analysis/R/02_compiled_results.rds")


############################################# Create Plots ####################################################
res <- readRDS(file = "Analysis/R/02_compiled_results.rds")

# Recode the priors using plain text (as characters)
res$priors <- recode(
  res$priors,
  ref = "Reference",
  bhm = "BHM",
  commensurate = "Commensurate",
  npp = "NPP",
  leap = "LEAP",
  `pp_0.2` = "PP_0.2",
  `pp_0.5` = "PP_0.5",
  `pp_0.8` = "PP_0.8",
  psipp = "PSIPP"
)
res$priors <- factor(
  res$priors,
  levels = c(
    "Reference", "BHM", "Commensurate",
    "PP_0.2", "PP_0.5", "PP_0.8",
    "LEAP", "NPP", "PSIPP"
  ),
  ordered = T
)

# Create a custom labeller function to handle the expression for priors
prior_labeller <- function(x) {
  sapply(x, function(val) {
    if (val == "PP_0.2") {
      return(TeX(sprintf(r'(PP ($a_{0} = %s$))', 0.2)))
    } else if (val == "PP_0.5") {
      return(TeX(sprintf(r'(PP ($a_{0} = %s$))', 0.5)))
    } else if (val == "PP_0.8") {
      return(TeX(sprintf(r'(PP ($a_{0} = %s$))', 0.8)))
    } else {
      return(val)
    }
  })
}


# Pivot longer for metrics
res.long <- res %>%
  pivot_longer(
    cols = c("DIC", "LOOIC"),
    names_to = "metric",
    values_to = "value"
  )

# Create the plot
p1 <- res.long %>%
  filter(J > 4) %>%
  ggplot(aes(x = J, y = value, linetype = metric)) +
  geom_point() +
  scale_x_continuous(breaks = seq(2, 15, by = 2)) +
  geom_line() +
  labs(x = 'J', linetype = 'Information Criteria') +
  facet_wrap(~ priors, labeller = labeller(priors = as_labeller(prior_labeller, default = label_parsed))) +
  labs(title = "DIC/LOOIC vs. Number of Intervals (J)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
p1

p2 = res %>%
  filter(J > 4) %>%
  ggplot(aes(x = J, y = logml)) +
  geom_point() +
  scale_x_continuous(breaks = seq(2, 15, by = 2))+
  geom_line() +
  labs(x = 'J', y = 'log(marginal likelihood)') +
  facet_wrap(~ priors, labeller = labeller(priors = as_labeller(prior_labeller, default = label_parsed))) +
  labs(title = "Log Marginal Likelihood v.s. Number of Intervals (J)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p2
