#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=
#SBATCH --job-name="kraken"
#SBATCH --output=kraken.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Modules loading
module purge # unloads all previous loads
module load miniconda3/default
source /mnt/big/miniconda3/condainit.sh
conda activate metawrap-env
module load metaWRAP/1.3
module load ncbi-blast/2.10.0+

# 3 Kraken module on both reads and the assembly
metawrap kraken2  -o KRAKEN \
                  -t 40 \
                  ASSEMBLY/final_assembly.fasta
                  CLEAN_READS/Elos*fastq \


# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

