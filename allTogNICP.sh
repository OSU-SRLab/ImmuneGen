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

#make a directory for the results
#mkdir /fs/ess/PAS0854/Active_projects/`echo $name`
#make a directory for the bams
mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/bams

#sbatch makeOrienmaf.sh `echo $reg` `echo $name` `echo $user`

#sbatch makeTCGAmaf.sh `echo $reg` `echo $name` `echo $user`

sbatch getCopy.sh `echo $name` `echo $user`

sbatch getFusion.sh `echo $name` `echo $user`

#load R
module load R/4.1.0-gnu9.1

echo "Making toSlice"

#make the list of files to slice
Rscript redoSlice.R `echo $name`

echo "Downloading"

#download the files
sbatch runDownloadEv.sh `echo $name` `echo $reg`

echo "Waiting: Download";
Rscript waiting.R;

cond=`squeue -u "$user" | grep 'runD' | wc -l`

#wait for first round of downloads to run
while [ $cond -gt 0 ];
do 
    echo "Waiting: Download";
    Rscript waiting.R;
    cond=`squeue -u "$user" | grep 'runD' | wc -l`;
done

#record the samples downloaded
cond=`ls /fs/ess/PAS0854/Active_projects/"$name"/bams | wc -l`
#record the total number of samples to be downloaded
lines=`cat /fs/ess/PAS0854/Active_projects/"$name"/toSlice.txt | wc -l`

#after the first round is run, check to see if there was a 500 error
while [ $cond -lt $lines ];
do
    #if there was a 500 error, remake the list to slice
    Rscript redoSlice.R `echo $name`;
    #retry the downloads
    sbatch runDownloadEv.sh `echo $name` `echo $reg`;
    echo "Waiting: Download";
    Rscript waiting.R;
    #wait for the next round of downloads to run
    cond=`squeue -u "$user" | grep 'runD' | wc -l`
    while [ $cond -gt 0 ];
    do 
            echo "Waiting: Download";
            Rscript waiting.R;
            cond=`squeue -u "$user" | grep 'runD' | wc -l`;
    done;
    #check to see if this round also encountered a 500 error
    cond=`ls /fs/ess/PAS0854/Active_projects/"$name"/bams | wc -l`;
    lines=`cat /fs/ess/PAS0854/Active_projects/"$name"/toSlice.txt | wc -l`;
done

#after everything is downloaded, rename the bams
Rscript renameBams.R `echo $name`

#rm /fs/ess/PAS0854/Active_projects/`echo $name`/CD40.PairedBams.txt

Rscript makePairedBams.R `echo $name`

module load samtools

#make a list of all the bams
ls /fs/ess/PAS0854/Active_projects/`echo $name`/bams/*.bam >> allBams.txt

#index all the bams
while IFS= read line;
do
    samtools index `echo $line`;
done < allBams.txt

#remove the list of the bams
rm allBams.txt

#make a directory for the delly results
mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/delly

#run delly on the bams
sbatch DellyAll.sh `echo $name`

#make gridss directory
mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/gridss

#run gridss
sbatch ./gridssFiles/GridssAll.sh `echo $reg` `echo $name` `echo $user`

#make a manta directory
mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/manta

#run manta
sbatch ./mantaFiles/MantaAll.sh `echo $reg` `echo $name` `echo $user`

cond=`squeue -u $user | grep 'All.sh' | wc -l`
while [ $cond -ge 0 ];
do
    #if we are, wait 30 min and check again
    echo "Waiting: SV callers";
    Rscript waiting.R;
    cond=`squeue -u $user | grep 'All.sh' | wc -l`;
done

#put everything into one data frame
sbatch test.sh `echo $reg` `echo $name` `echo $user`
sbatch remove.sh
