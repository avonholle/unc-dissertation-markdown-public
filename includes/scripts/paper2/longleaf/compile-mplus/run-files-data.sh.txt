#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -o  run-mplus-prep.Rout
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=vonholle@email.unc.edu # send-to address

module add r
module add rstudio
R CMD BATCH run-mplus-prep.R