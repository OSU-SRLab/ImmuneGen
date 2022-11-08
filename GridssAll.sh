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

#go line by line through the paired bam file
while IFS= read line;
do
    cond=`squeue -u $user | wc -l`
    #check that we are not going to exceed 1000 jobs
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Gridss";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done;
    #run gridss
    A=`echo $line | cut -d ";" -f1`; 
    B=`echo $line | cut -d ";" -f2`;
    sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/gridssFiles/runGridss1ata.sh `echo $A` `echo $B` `echo $name`;
done < /fs/ess/PAS0854/Active_projects/`echo $name`/CD40.PairedBams.txt

while IFS= read line;
do
    cond=`squeue -u $user | wc -l`
    #check that we are not going to exceed 1000 jobs
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Gridss";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done;
    #run gridss
    A=`echo $line | cut -d ";" -f1`; 
    B=`echo $line | cut -d ";" -f2`;
    sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/gridssFiles/runGridss1ataOrien.sh `echo $A` `echo $B` `echo $name`;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/AllPairedBamsOrien.txt

while IFS= read line;
do
    cond=`squeue -u $user | wc -l`
    #check that we are not going to exceed 1000 jobs
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Gridss";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done;
    #run gridss
    A=`echo $line | cut -d ";" -f1`; 
    B=`echo $line | cut -d ";" -f2`;
    sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/gridssFiles/runGridss1ataOrien.sh `echo $A` `echo $B` `echo $name`;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/CD40.PairedBamsOrien2.txt
