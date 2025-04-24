## Set up additional simulation scenarios (nevents = 600)
library(dplyr)

## Directory to save simulation scenarios
save.dir <- 'Sims_supp/R'

## Simulation scenario
nevents    <- 600
n0events   <- 150
cens.prop  <- 0.4 ## censoring proportion
case       <- c("UNEXCH", "BHM", "PP")
grid       <- expand.grid(
  nevents = nevents, n0events = n0events,
  cens.prop = cens.prop, case = case
  , stringsAsFactors = FALSE
)

ndatasets <- 10000
each.cl   <- 10 ## how many data sets to run on a single node
ncl       <- ceiling(ndatasets / each.cl) ## how many nodes per data set

## Repeat each row of the grid ncl times
grid <- grid[rep(1:nrow(grid), each = ncl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1) ## 3000 rows
saveRDS(grid, file = file.path(save.dir, "grid.rds"))

set.seed(1)
# generate seeds for each data set
seeds_list <- sample(seq_len(1000 * ndatasets), ndatasets, replace = FALSE)
saveRDS(seeds_list, file = file.path(save.dir, "seeds_list.rds"))
