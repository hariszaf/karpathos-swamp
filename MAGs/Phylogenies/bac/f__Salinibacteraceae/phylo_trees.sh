#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --nodelist=
#SBATCH --ntasks-per-node=10
#SBATCH --mem=55000
#SBATCH --job-name="f__Salinibacteraceae"
#SBATCH --output=iqtree.output
#SBATCH --mail-user=haris.zafr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Input:
# 1. wget from gtdb all the entries of the taxonomy group of interest (all_members.tsv)
# 2. a file with the MAGs names of interest ()


module purge

# Get accession numbers of the representative genomes
awk -F"\t" '$5=="True" {print $1 }' all_members.tsv  > repr_accessions_numbers.tsv

# From the GTDB MSA get the alignments of the representative genomes
while IFS= read -r line; 
do
    grep -A 1 $line /mnt/big/miniconda3/envs/protologger/share/gtdbtk-1.5.0/db/pplacer/gtdb_r202_bac120.refpkg/bac120_msa_reps_r202.faa >> repr_seqs.aln ; 
done < repr_accessions_numbers.tsv 

# From the GTDB-Tk alignment of our MAGs, get the alignments of the novel taxa of interest
while IFS= read -r line; 
do
    grep -A 1 -w $line /home1/christina/Elos_meta/Curtain3/metaWRAP/GTDB_TK_bins_classification/elos_70_5.bac120.user_msa.fasta  >> users_msa.aln ;
done < my_mags.tsv

# Contatenate the 2 alignment sections
cat repr_seqs.aln users_msa.aln  > all.aln

# Make a tree out of them
/mnt/big/Phylogeny/iqtree-2.0-rc1-Linux/bin/iqtree  -s all.aln \
                                                    -m MFP \
                                                    --prefix investigation_bs \
                                                    -alrt 1000 \
                                                    -B 1000 \
                                                    -T AUTO



