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

#go line by line through the paired bam file
while IFS= read line;
do
    #check that we are not going to exceed 1000 jobs
    cond=`squeue -u $user | wc -l`
    while [ $cond -ge 999 ];
    do
        #if we are, wait 30 min and check again
        echo "Waiting: Manta";
        Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
        cond=`squeue -u $user | wc -l`;
    done; 
    #run manta
    A=`echo $line | cut -d ";" -f1`;
    B=`echo $line | cut -d ";" -f2`;
    C=`echo $line | cut -d ";" -f3`;
    mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
    /fs/ess/PAS0854/Software/manta-1.4.0.centos6_x86_64/bin/configManta.py \
    --normalBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $C` \
    --tumorBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $B` \
    --referenceFasta /fs/ess/PAS0854/Reference_Data/GRCh38.d1.vd1.fa \
    --runDir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
    sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/mantaFiles/runManta1ata.sh `echo $name` `echo $A`;
done < /fs/ess/PAS0854/Active_projects/`echo $name`/CD40.PairedBams.txt

# #go line by line through the paired bam file
# while IFS= read line;
# do
#     #check that we are not going to exceed 1000 jobs
#     cond=`squeue -u $user | wc -l`
#     while [ $cond -ge 999 ];
#     do
#         #if we are, wait 30 min and check again
#         echo "Waiting: Manta";
#         Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
#         cond=`squeue -u $user | wc -l`;
#     done; 
#     #run manta
#     A=`echo $line | cut -d ";" -f1`;
#     B=`echo $line | cut -d ";" -f2`;
#     C=`echo $line | cut -d ";" -f3`;
#     mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
#     /fs/ess/PAS0854/Software/manta-1.4.0.centos6_x86_64/bin/configManta.py \
#     --normalBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $C` \
#     --tumorBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $B` \
#     --referenceFasta /fs/ess/PAS0854/Reference_Data/ORIEN/hs38DH.fa \
#     --runDir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
#     sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/mantaFiles/runManta1ata.sh `echo $name` `echo $A`;
# done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/CD40.PairedBamsOrien2.txt

# #go line by line through the paired bam file
# while IFS= read line;
# do
#     #check that we are not going to exceed 1000 jobs
#     cond=`squeue -u $user | wc -l`
#     while [ $cond -ge 999 ];
#     do
#         #if we are, wait 30 min and check again
#         echo "Waiting: Manta";
#         Rscript /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/waiting.R;
#         cond=`squeue -u $user | wc -l`;
#     done; 
#     #run manta
#     A=`echo $line | cut -d ";" -f1`;
#     B=`echo $line | cut -d ";" -f2`;
#     C=`echo $line | cut -d ";" -f3`;
#     mkdir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
#     /fs/ess/PAS0854/Software/manta-1.4.0.centos6_x86_64/bin/configManta.py \
#     --normalBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $C` \
#     --tumorBam /fs/ess/PAS0854/Active_projects/`echo $name`/bams/`echo $B` \
#     --referenceFasta /fs/ess/PAS0854/Reference_Data/ORIEN/hs38DH.fa \
#     --runDir /fs/ess/PAS0854/Active_projects/`echo $name`/manta/`echo $A`.manta;
#     sbatch /fs/project/PAS0854/Active_projects/SV_SNV_CompBatch/mantaFiles/runManta1ata.sh `echo $name` `echo $A`;
# done < /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/lists/AllPairedBamsOrien.txt
