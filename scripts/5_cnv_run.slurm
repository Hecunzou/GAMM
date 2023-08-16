#!/bin/bash
#
module load apps/jags/4.3.1
module load apps/R/4.1.0
echo ${SLURM_JOB_NODELIST}
pwd
Rscript ./cnv_smart.R

echo $?
echo end on $(date)
##sbatch cnv_run.slurm

