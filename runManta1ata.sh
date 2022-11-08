#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=12:00:00
#SBATCH --nodes=1

module load python/2.7-conda5.2

python /fs/ess/PAS0854/Active_projects/`echo $1`/manta/$2.manta/runWorkflow.py -m local;