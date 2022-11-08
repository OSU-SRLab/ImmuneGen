#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

module load R/4.1.0-gnu9.1
Rscript download_TCGAEve.R `echo $1` `echo $2` > download_TCGAEve.Rout
