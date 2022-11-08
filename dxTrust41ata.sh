#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=72:00:00
#SBATCH --nodes=1

#take in arguments from user
#name is the name of the gene
name=`echo $2`

cd /fs/ess/PAS0854/Active_projects/`echo $name`/trust4

#token=$(<'/fs/ess/PAS0854/Reference_Data/GDC_Token/gdc-user-token.2022.09.08.txt')

#curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/data/'`echo $1` -o `echo $1`.bam;

run-trust4 --abnormalUnmapFlag -b ./`echo $1`.bam -f /fs/ess/PAS0854/Software/TRUST4/hg38_bcrtcr.fa --ref /fs/ess/PAS0854/Software/TRUST4/human_IMGT+C.fa
