#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="anvSum"
#SBATCH --output=summary.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

anvi-summarize -p SAMPLES-MERGED-1/PROFILE.db -c swamp-assembly.db -C METWR_REFINED_BIN_COLLECTION -o FIRST_SUMMARY

module purge

