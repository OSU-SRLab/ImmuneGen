#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=10:00:00
#SBATCH --nodes=1

while IFS= read line;
do
sbatch dxTrust4Orien1ata.sh $line;
done < /fs/project/PAS0854/Raven/Immune_Checkpoints/CD40.RNAseqBams.Orien1.txt
