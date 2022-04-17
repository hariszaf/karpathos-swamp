#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=40
#SBATCH --mem=
#SBATCH --job-name="refinement"
#SBATCH --output=bin_refinement.output
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


# 5. CONSOLIDATE BINS
metawrap bin_refinement -o BIN_REFINEMENT \
			-t 40 \
			-A INITIAL_BINNING/concoct_bins/ \
			-B INITIAL_BINNING/maxbin2_bins/ \
			-c 70 \
			-x 5

## 6. BLOBOLOGY MODULE						[Currently not working due to blastn issue]
#metawrap blobology -a ASSEMBLY/final_assembly.fasta \
#		   -t 20 \
#		   -o BLOBOLOGY \
#		   --bins BIN_REFINEMENT/metawrap_70_5_bins \
#		   CLEAN_READS/Elos*fastq



# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

