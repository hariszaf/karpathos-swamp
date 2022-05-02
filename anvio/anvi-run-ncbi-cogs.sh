#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="anv-cogs"
#SBATCH --output=met_estimate_internal_genomes.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1


anvi-run-ncbi-cogs -c ../swamp-assembly.db  --num-threads 40 --temporary-dir-path /home1/christina/ANVIO_TMP

module purge

