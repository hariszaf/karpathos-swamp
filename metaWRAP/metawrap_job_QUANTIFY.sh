#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=40
#SBATCH --mem=
#SBATCH --job-name="quantify"
#SBATCH --output=quantify.output
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

# 7. ABUNDANCE OF GENOMES ACCROSS SAMPLES
metawrap quant_bins -b BIN_REFINEMENT/metawrap_70_5_bins \
		    -o QUANT_BINS \
		    -t 40 \
		    -a ASSEMBLY/final_assembly.fasta \
		    CLEAN_READS/Elos*fastq


# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

