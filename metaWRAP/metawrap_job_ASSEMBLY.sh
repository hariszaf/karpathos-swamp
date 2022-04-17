#!/bin/bash -l

#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=20
#SBATCH --mem=
#SBATCH --job-name="coassembly"
#SBATCH --output=assembly_all.output
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


## 2. ASSEMBLY
# The -m parameter denotes the RAM GB to allocate
metawrap assembly --megahit \
                  -m 1000 \
                  -t 20 \
		  -1 CLEAN_READS/ALL_READS_1.fastq \
		  -2 CLEAN_READS/ALL_READS_2.fastq \
		  -o ASSEMBLY


# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

