#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=40
#SBATCH --mem=
#SBATCH --job-name="metabat"
#SBATCH --output=bin_metabat2.output
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


# 4. BIN THE (CO-)ASSEMBLY
# 4.c METABAT2
metawrap binning -o INITIAL_BINNING \
		 -t 40 \
		 -a ASSEMBLY/final_assembly.fasta \
		 --metabat2 CLEAN_READS/Elos*fastq \
		 --run-checkm

# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

