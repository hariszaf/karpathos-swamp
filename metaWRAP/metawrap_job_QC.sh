#!/bin/bash -l

#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=20
#SBATCH --mem=
#SBATCH --job-name="QC"
#SBATCH --output=quality_control.output
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

## 1. READQC
for F in /home1/christina/Elos_meta/CONCATENATED_RAW_DATA/*_1.fastq; do 
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap read_qc -1 $F -2 $R -t 20 -o READ_QC/$SAMPLE
done




# Modules unload
module unload miniconda3/default
module unload metaWRAP/1.3

