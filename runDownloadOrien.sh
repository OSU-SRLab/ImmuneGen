#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

#take in arguments from user
#reg is the region to slice
reg=`echo $1`
#name is the name of the gene
name=`echo $2`
#user is the username
user=`echo $3`

module load R/4.0.2-gnu9.1

mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/empty

cd /fs/ess/PAS0854/Active_projects/`echo $name`/empty

sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/orienFiles/makeBatchesDxOrien.sh `echo $reg` `echo $name`

echo "Waiting: OrienDownload"
Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R

cond=`squeue -u "$user" | grep 'Dx' | wc -l`
while [ $cond -gt 0 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Orien";
        Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | grep 'Dx' | wc -l`
    done;
