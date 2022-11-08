#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

name=`echo $1`
user=`echo $2`

mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/copy

module load R/4.1.0-gnu9.1

while IFS= read line;
do
    cond=`squeue -u $user | wc -l`
    #check that we are not going to exceed 1000 jobs
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Copy";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done;
    sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/getCopy1ata.sh `echo $name` `echo $line`;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/copy_number.txt

echo "Waiting: Copy Download"
Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R

cond=`squeue -u "$user" | grep 'getC' | wc -l`
while [ $cond -gt 1 ];
do
    #if we are, wait 30 min and check again
    echo "Waiting: Copy Download";
    Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
    cond=`squeue -u $user | grep 'getC' | wc -l`;
done

echo "Waiting: Making copy.num.txt"

cd /fs/ess/PAS0854/Active_projects/`echo $name`/copy

ls > ls.txt

echo -e -n "Sample\t" >> /fs/ess/PAS0854/Active_projects/`echo $name`/copy.num.txt
head -1 00031471-9649-4d9e-b91b-cb5139e41883.tsv >> /fs/ess/PAS0854/Active_projects/`echo $name`/copy.num.txt

while IFS= read line; 
do 
    A=`echo $line`; 
    tail -n+2 /fs/ess/PAS0854/Active_projects/`echo $name`/copy/`echo $A` | grep `echo $name` | awk -v prefix="$A" '{print prefix "\t" $0}'>> /fs/ess/PAS0854/Active_projects/`echo $name`/copy.num.txt;
done < /fs/ess/PAS0854/Active_projects/`echo $name`/copy/ls.txt

rm -r /fs/ess/PAS0854/Active_projects/`echo $name`/copy
