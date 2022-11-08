#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=10:00:00
#SBATCH --nodes=1

#take in arguments from user
name=`echo $1`
user=`echo $2`

#mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/trust4

cd /fs/ess/PAS0854/Active_projects/`echo $name`

module load R/4.1.0-gnu9.1

#Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/makeRNAlist.R

cond=`squeue -u $user | wc -l`
while IFS= read line;
do
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Download";
        Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`
    done;
sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/dxTrust41ata.sh `echo $line` `echo $name`;
cond=`squeue -u $user | wc -l`
done < RNAseqBams.TCGA.txt

cond=`squeue -u $user | wc -l`
while IFS= read line;
do
while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Download";
        Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`
    done;
sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/dxTrust4Orien1ata.sh `echo $line` `echo $name`;
cond=`squeue -u $user | wc -l`
done < RNAseqBams.Orien.txt

Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R
cond=`squeue -u $user | grep 'dxTr' | wc -l`

while [ $cond -gt 0 ];
do
    echo "Waiting: Download";
    Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R; 
    cond=`squeue -u $user | grep 'dxTr' | wc -l`;
done

rm slurm-*

Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/Trust4Analysis.R `echo $name`
