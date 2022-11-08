#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=72:00:00
#SBATCH --nodes=1

module load python
module load samtools

while IFS= read line; 
do A=`echo $line | cut -d ";" -f1`; 
B=`echo $line | cut -d ";" -f2`;
C=`echo $line | cut -d ";" -f3`;
delly call -o /fs/ess/PAS0854/Active_projects/`echo $1`/delly/`echo $A`.bcf -g /fs/ess/PAS0854/Reference_Data/GRCh38.d1.vd1.fa /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $B` /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $C`;
done < /fs/ess/PAS0854/Active_projects/`echo $1`/CD40.PairedBams.txt

# while IFS= read line; 
# do A=`echo $line | cut -d ";" -f1`; 
# B=`echo $line | cut -d ";" -f2`;
# C=`echo $line | cut -d ";" -f3`;
# samtools index /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $B`;
# samtools index /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $C`;
# delly call -o /fs/ess/PAS0854/Active_projects/`echo $1`/delly/`echo $A`.bcf -g /fs/ess/PAS0854/Reference_Data/ORIEN/hs38DH.fa /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $B` /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $C`;
# done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/AllPairedBamsOrien.txt

# while IFS= read line; 
# do A=`echo $line | cut -d ";" -f1`; 
# B=`echo $line | cut -d ";" -f2`;
# C=`echo $line | cut -d ";" -f3`;
# samtools index /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $B`;
# samtools index /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $C`;
# delly call -o /fs/ess/PAS0854/Active_projects/`echo $1`/delly/`echo $A`.bcf -g /fs/ess/PAS0854/Reference_Data/ORIEN/hs38DH.fa /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $B` /fs/ess/PAS0854/Active_projects/`echo $1`/bams/`echo $C`;
# done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/CD40.PairedBamsOrien2.txt
