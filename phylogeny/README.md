



To get the single-copy genes for each MAG from the anvio' database: 


```
while IFS= read -r gene; do  anvi-get-sequences-for-hmm-hits -c swamp-assembly.db -p PROFILE.db -C METWR_REFINED_BIN_COLLECTION -o best_single_copy_aa_concat_$gene.fa --hmm-source Bacteria_71 --gene-names $gene --return-best-hit --get-aa-sequences  --concatenate  ; done < common-genes
```

