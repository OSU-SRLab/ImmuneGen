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

#Now go to the orien files and read them in to Orien.SNV.maf
cd /fs/ess/PAS0854/Active_projects/TCGA_novel_icktps/OrienMafs
A=`head -1 A35748_st_t_markdup_recalibrated_Haplotyper.ft.M2GEN.PoN.v2.vcf.out.hg38_multianno.txt`
echo -e -n "Sample\t" >> /fs/ess/PAS0854/Active_projects/`echo $name`/Orien.SNV.maf
echo "$A" >> /fs/ess/PAS0854/Active_projects/`echo $name`/Orien.SNV.maf
while IFS= read line; 
do
cat ./"$line" | grep `echo $name` | awk -v prefix="$line" '{print prefix "\t" $0}' >> /fs/ess/PAS0854/Active_projects/`echo $name`/Orien.SNV.maf;
done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/OrienMafs.txt
