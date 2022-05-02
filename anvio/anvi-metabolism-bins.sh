#!/bin/bash -l

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="metabolism"
#SBATCH --output=met_estimate_internal_genomes.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

# This provides the overview file of what modules are in each MAG; required step 
# to run the internal-genomes
#anvi-estimate-metabolism -c swamp-assembly.db -p SAMPLES-MERGED-1/PROFILE.db -C MAGs

# Here is the per bin case
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format



module purge
