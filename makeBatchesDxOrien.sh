#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=10:00:00
#SBATCH --nodes=1

#put the lists into a master list called use.txt
while IFS= read line; do
    sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/orienFiles/DxSliceOrien.sh $line `echo $1` `echo $2`;
done < /fs/ess/PAS0854/Active_projects/TCGA_novel_icktps/toDownload/use.txt
#run DxSliceOrien on each list
