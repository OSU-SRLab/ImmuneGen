#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=1

name=`echo $1`
UUID=`echo $2`

token=$(</fs/ess/PAS0854/Reference_Data/GDC_Token/gdc-user-token.2022.10.19.txt)

curl -o "/fs/ess/PAS0854/Active_projects/${name}/fusion/${UUID}.tsv" --remote-header-name --header "X-Auth-Token: $token" "https://api.gdc.cancer.gov/data/${UUID}"