#!/bin/bash -l

#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="singProf"
#SBATCH --output=single_profile_12.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

anvi-profile -i MAPPING/Elos1.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos1COGs
#anvi-profile -i MAPPING/Elos2.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos2
#anvi-profile -i MAPPING/Elos3.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos3
#anvi-profile -i MAPPING/Elos7.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos7
#anvi-profile -i MAPPING/Elos10.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos10
#anvi-profile -i MAPPING/Elos12.bam -c swamp-assembly.db --num-threads 40 --sample-name profileElos12

module purge

