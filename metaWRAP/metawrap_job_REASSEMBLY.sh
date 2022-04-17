#!/bin/bash -l

#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="reassembly"
#SBATCH --output=reassembly.output
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
module load SPAdes/3.14.0

# 8. RE-ASSEBMLE
# REMEMBER: Rename Elos-1-L7 to ALL_READS_ when you will be working with multiple samples
metawrap reassemble_bins -o BIN_REASSEMBLY \
			 -1 CLEAN_READS/ALL_READS_1.fastq \
			 -2 CLEAN_READS/ALL_READS_2.fastq \
			 -b BIN_REFINEMENT/metawrap_70_5_bins/ \
			 -t 20 \
			 -m 1500 \
			 -c 70 \
			 -x 5 \

# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

