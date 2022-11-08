#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

name=`echo $1`
user=`echo $2`

mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/fusion

module load R/4.1.0-gnu9.1

while IFS= read line;
do
    cond=`squeue -u $user | wc -l`
    #check that we are not going to exceed 1000 jobs
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Fusion";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done;
    sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/getFusion1ata.sh `echo $name` `echo $line`;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/fusion.txt

echo "Waiting: Fusion Download"
Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R

cond=`squeue -u "$user" | grep 'getF' | wc -l`
while [ $cond -gt 1 ];
do
    #if we are, wait 30 min and check again
    echo "Waiting: Fusion Download";
    Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
    cond=`squeue -u $user | grep 'getF' | wc -l`;
done

echo "Waiting: Making Fusion.txt"

cd /fs/ess/PAS0854/Active_projects/`echo $name`/fusion

ls > ls.txt

echo -e -n "Sample\t" >> /fs/ess/PAS0854/Active_projects/`echo $name`/fusion.txt

head -1 000377d3-6a66-4552-9c4b-5cd4487adc9a.tsv >> /fs/ess/PAS0854/Active_projects/`echo $name`/fusion.txt

while IFS= read line; 
do 
    A=`echo $line`; 
    tail -n+2 /fs/ess/PAS0854/Active_projects/`echo $name`/fusion/`echo $A` | grep `echo $name` | awk -v prefix="$A" '{print prefix "\t" $0}'>> /fs/ess/PAS0854/Active_projects/`echo $name`/fusion.txt;
done < /fs/ess/PAS0854/Active_projects/`echo $name`/fusion/ls.txt

rm -r /fs/ess/PAS0854/Active_projects/`echo $name`/fusion
