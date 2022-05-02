#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="hmm-db"
#SBATCH --output=hmms-db.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

anvi-run-hmms -T 20 -c swamp-assembly.db

module purge

