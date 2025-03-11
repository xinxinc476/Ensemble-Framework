#!/bin/bash

#SBATCH --job-name=Analysis
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000m
#SBATCH	--array=[1-196]
#SBATCH --output=/work/users/x/i/xinxinc/super_prior/Analysis/Log/slurmLogFiles%a.out
#SBATCH --error=/work/users/x/i/xinxinc/super_prior/Analysis/Error/%a.err

## add R module
module add gcc/11.2.0
module add r/4.3.1

R CMD BATCH --no-restore ~/projects/super_prior/Analysis/R/02_fit_PWE_CurePWE.R /work/users/x/i/xinxinc/super_prior/Analysis/Rout/$SLURM_ARRAY_TASK_ID.Rout
