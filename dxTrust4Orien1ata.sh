#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=72:00:00
#SBATCH --nodes=1

#take in arguments from user
#name is the name of the gene
name=`echo $2`

cd /fs/ess/PAS0854/Active_project/`echo $name`/trust4

module load samtools

#dx download $1

A=`echo $1 | rev | cut -d "/" -f1 | rev`
B=`echo $A | rev | cut -d "." -f2- | rev`

samtools view -b -T /fs/project/PAS0854/Reference_Data/ORIEN/hs38DH.fa -o `echo $B`.bam `echo $A`

rm `echo $A`

run-trust4 --abnormalUnmapFlag -b ./`echo $B`.bam -f /fs/ess/PAS0854/Software/TRUST4/hg38_bcrtcr.fa --ref /fs/ess/PAS0854/Software/TRUST4/human_IMGT+C.fa
