#!/bin/bash -l

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="kofam"
#SBATCH --output=kofam.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

export TMPDIR="/home1/christina/ANVIO_TMP"
anvi-run-kegg-kofams -c swamp-assembly.db  -T 40

module purge
