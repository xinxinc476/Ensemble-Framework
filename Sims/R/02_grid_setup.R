## Set up simulation scenarios for simulation
remove(list = ls())

library(dplyr)

## Directory to save simulation scenarios
save.dir <- 'Sims/R'

## Simulation scenario
n    <- seq(150, 500, by = 50)
#n0   <- 250
case <- c("non-exchangeable", "PP", "BHM", "LEAP", "PSIPP")
grid <- expand.grid(
  n = n, case = case
  , stringsAsFactors = FALSE
)

ndatasets <- 10000  ## total number of data sets
each.cl   <- 400 ## how many data sets to run on a single node
ncl       <- ceiling(ndatasets / each.cl) ## how many nodes per data set

## Repeat each row of the grid ncl times (want nrow(grid) <= 800 ideally)
grid <- grid[rep(1:nrow(grid), each = ncl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1)
saveRDS(grid, file = file.path(save.dir, "grid.rds"))

set.seed(1)
# generate seeds for each data set
seeds_list <- sample(seq_len(1000 * ndatasets), ndatasets, replace = FALSE)
saveRDS(seeds_list, file = file.path(save.dir, "seeds_list.rds"))
