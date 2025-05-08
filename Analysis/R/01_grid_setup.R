## Set up grid for real data analysis
library(dplyr)

## Directory to save simulation scenarios
save.dir <- 'Analysis/R'

## Obtain scenarios for computation
model <- c("PWE", "CurePWE")
prior <- c('ref', 'pp', 'psipp', 'bhm', 'cp', 'leap', 'npp')
J      <- 2:9 # number of intervals for the PWE model
# 9 is the largest number such that there is at least one event within each interval
# for the stratified current and historical data sets
arm    <- c(0, 1) # we will stratify the current data by treatment arm for analysis
grid   <- expand.grid(
  J = J, arm = arm, prior = prior, model = model
  , stringsAsFactors = FALSE
) ## nrow = 224

saveRDS(grid, file = file.path(save.dir, "grid.rds"))
