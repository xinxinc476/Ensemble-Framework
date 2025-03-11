## Set up grid for real data analysis
library(dplyr)

## Directory to save simulation scenarios
save.dir <- 'Analysis/R'

## Obtain scenarios for computation
models <- c("PWE", "CurePWE")
priors <- c('ref', 'pp', 'psipp', 'bhm', 'cp', 'leap', 'npp')
J      <- 2:15 # number of intervals for the PWE model
grid   <- expand.grid(
  J = J, priors = priors, models = models
  , stringsAsFactors = FALSE
) ## nrow = 196

saveRDS(grid, file = file.path(save.dir, "grid.rds"))
