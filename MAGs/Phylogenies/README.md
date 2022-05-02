In each folder you will find a `all_members.tsv` file
including all the GTDB genomes that are part of the taxon under study. 

A `my_mags.tsv` file will be also present in every folder, including the MAGs 
of this study that were assigned within the taxon under study.

Using the `phylo_trees.sh` script the phylogenetic tree for each case was built 
getting the corresponding line from the GTDB msa file for the genomes already 
on GTDB and those of our MAGs from their assignment using GTDB-Tk. 

Finally, the `get_representatives.py` script was used to calculate the POCP values. 



