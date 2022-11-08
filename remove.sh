#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

cd /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch
rm slurm*

cd /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/gridssFiles

rm -r *.gridss.working
rm gridss.*
rm slurm-*

cd /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/mantaFiles

rm slurm-*
