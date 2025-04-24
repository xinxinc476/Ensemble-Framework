## Set up simulation scenarios
library(dplyr)

## Directory to save simulation scenarios
save.dir <- 'Sims/R/case_Exch'

## Simulation scenario
nevents    <- c(50, 150, 300)
n0events   <- 150
cens.prop  <- c(0.2, 0.4, 0.6) ## censoring proportion
case       <- "PP"
grid       <- expand.grid(
  nevents = nevents, n0events = n0events,
  cens.prop = cens.prop, case = case
  , stringsAsFactors = FALSE
)

ndatasets <- 10000
each.cl   <- 25 ## how many data sets to run on a single node
ncl       <- ceiling(ndatasets / each.cl) ## how many nodes per data set

## Repeat each row of the grid ncl times (want nrow(grid) <= 800 ideally)
grid <- grid[rep(1:nrow(grid), each = ncl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1) ## 3600 rows
saveRDS(grid, file = file.path(save.dir, "grid.rds"))

set.seed(1)
# generate seeds for each data set
seeds_list <- sample(seq_len(1000 * ndatasets), ndatasets, replace = FALSE)
saveRDS(seeds_list, file = file.path(save.dir, "seeds_list.rds"))
