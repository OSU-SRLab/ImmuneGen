#!/bin/sh
#SBATCH --account=PAS0854
#SBATCH --time=96:00:00
#SBATCH --nodes=5

#take in arguments from user
#reg is the region to slice
reg=`echo $1`
#name is the name of the gene
name=`echo $2`
#user is the username
user=`echo $3`

module load R/4.1.0-gnu9.1

# Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/makeMaster.R `echo $name`

# Rscript /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/makeControls.R `echo $name`

# sbatch /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/dxTrust4.sh `echo $name` `echo $user`

cd /fs/ess/PAS0854/Active_projects/SV_SNV_CompBatch/makeGraphs

Rscript assignFunction.R `echo $name`

# Rscript makeBreakpoints.R `echo $name` `echo $reg`

# Rscript makeSNVmap.R `echo $name` `echo $reg`

# Rscript makeCancerType.R `echo $name`

# Rscript makeExpression.R `echo $name`
