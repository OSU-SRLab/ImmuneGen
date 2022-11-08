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

#now get the SNV data for the gene
#go to the TCGA files first
cd /fs/ess/PAS0854/Active_projects/TCGA_novel_icktps/TCGA.SomaticMut.Mafs

#read them into a file that is called TCGA.SNV.maf
zcat 00020280-6040-4071-b30e-b882a3dfbc2b.wxs.Pindel.aliquot.maf.gz | grep -v '#' | head -1 >> /fs/ess/PAS0854/Active_projects/`echo $name`/TCGA.SNV.maf
while IFS= read line; 
do A=`echo $line | cut -d "," -f4 | tr -d '"'`;
B=`echo $line | cut -d "," -f2 | tr -d '"'`; 
zcat ./`echo $B` | grep -v '#' | awk 'NR>1' | grep `echo $name` >> /fs/ess/PAS0854/Active_projects/`echo $name`/TCGA.SNV.maf;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/Maf.Manifest.txt