#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="builddb"
#SBATCH --output=db-building.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

anvi-gen-contigs-database -f contigs-fixed.fa \
                          --num-threads 20 \
                          --project-name swamp-assembly \
                          -o swamp-assembly.db

module purge

