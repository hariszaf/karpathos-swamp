#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=
#do not change the above
#change the following properly
#SBATCH --job-name="seqshmms"
#SBATCH --output=get_seqs_hmms.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

module purge
module load anvio/7.1

/mnt/big/Metagenomics/anvio_v6.1/bin/anvi-get-sequences-for-hmm-hits -c swamp-assembly.db \
                                --hmm-source Bacteria_71 \
                                -o genes-fasta-bac

module purge

