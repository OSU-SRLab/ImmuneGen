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

#sbatch DellyAll.sh `echo $name`

sbatch ./gridssFiles/GridssAll.sh `echo $reg` `echo $name` `echo $user`

sbatch ./mantaFiles/MantaAll.sh `echo $reg` `echo $name` `echo $user`

cond=`squeue -u $user | grep 'All.sh' | wc -l`
while [ $cond -ge 0 ];
do
    #if we are, wait 30 min and check again
    echo "Waiting: SV callers";
    Rscript waiting.R;
    cond=`squeue -u $user | grep 'All.sh' | wc -l`
done;

#put everything into one data frame
sbatch test.sh `echo $reg` `echo $name` `echo $user`
sbatch remove.sh