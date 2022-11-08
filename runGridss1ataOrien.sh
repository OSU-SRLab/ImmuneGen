#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=3:00:00
#SBATCH --ntasks=16

singularity exec /fs/ess/PAS0854/Eric/gridss_latest.sif gridss /fs/ess/PAS0854/Active_projects/`echo $3`/bams/`echo $2` \
-r /fs/ess/PAS0854/Reference_Data/ORIEN/hs38DH.fa \
-t 8 \
-o /fs/ess/PAS0854/Active_projects/`echo $3`/gridss/`echo $1`.gridss.vcf.gz