#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

module load samtools

#read the list $1 into the while loop

while IFS= read line; do
#download cram
    dx download `echo $line`.cram;
#index cram
    samtools index `echo $line`.cram;
#slice cram
    samtools view -h -o /fs/ess/PAS0854/Active_projects/`echo $3`/bams/`echo $line`_sliced.cram `echo $line`.cram `echo $2`;
#remove the unsliced cram and the index
    rm `echo $line`.cram;
    rm `echo $line`.cram.crai;
done < /fs/ess/PAS0854/Active_projects/TCGA_novel_icktps/toDownload/$1
