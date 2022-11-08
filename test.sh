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

module load bcftools

cd /fs/ess/PAS0854/Active_projects/`echo $name`/manta

ls -d *.manta >> mantalist.txt
echo -e -n "Sample\t" >> ../MantaAllTog.txt
bcftools view ./TCGA-55-8620.manta/results/variants/candidateSV.vcf.gz | grep '#' | tail -1 >> ../MantaAllTog.txt
while IFS= read line; 
do A=`echo $line`; 
bcftools view ./`echo $A`/results/variants/candidateSV.vcf.gz | grep -v '#' | awk -v prefix="$A" '{print prefix "\t" $0}'>> ../MantaAllTog.txt;
done < mantalist.txt

rm mantalist.txt

#go to each file with that name, read the results out, and put it into a file

cd /fs/ess/PAS0854/Active_projects/`echo $name`/delly

echo -e -n "Sample\t" >> ../DellyAllTog.txt
bcftools view ./01A-10A-TCGA-02-2483.bcf | grep '#' | tail -1 >> ../DellyAllTog.txt
while IFS= read line; do A=`echo $line`; bcftools view `echo $A`.bcf | grep -v '#' | awk -v prefix="$A" '{print prefix "\t" $0}'>> ../DellyAllTog.txt; done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/dellylist.txt

cd /fs/ess/PAS0854/Active_projects/`echo $name`/gridss

echo -e -n "Sample\t" >> ../GridssAllTog.txt
bcftools view ./01A-10A-TCGA-02-2483.gridss.vcf.gz | grep '#' | tail -1 >> ../GridssAllTog.txt
while IFS= read line; do A=`echo $line`; bcftools view `echo $A`.gridss.vcf.gz | grep -v '#' | awk -v prefix="$A" '{print prefix "\t" $0}'>> ../GridssAllTog.txt; done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/dellylist.txt

module load R/4.1.0-gnu9.1
cd /fs/ess/PAS0854/Active_projects/`echo $name`
Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/make1file.R `echo $name`

sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/makeGraphs.sh `echo $reg` `echo $name` `echo $user`
