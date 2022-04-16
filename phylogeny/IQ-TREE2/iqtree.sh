#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --mem=
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="iqtree"
#SBATCH --output=findModel.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --requeue

# Get best model
#/mnt/big/Phylogeny/iqtree-2.0-rc1-Linux/bin/iqtree -s final_aa_trimmed_per_gene.aln -m MFP 

/mnt/big/Phylogeny/iqtree-2.0-rc1-Linux/bin/iqtree -s final_aa_trimmed_per_gene.aln \
                                                   -m LG+R10 \
                                                   -alrt 1000 \
                                                   -B 1000 \
                                                   -T AUTO




