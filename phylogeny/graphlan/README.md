`--pad` pad_in          the distance between the most external graphical
                        element and the border of the image


To get the clade names and annotate the tree with the bootstrap values, we ran: 

```
more mags_phylogeny.xml | grep "<name>" | grep -v "bin" | sed 's/\s//g' | sed 's/<name>// ; s/<\/name>//' > bs_clade_names.txt
```


To build a `graphlan` figure, you need to run: 

```bash
./build_annotion_file.py
graphlan_annotate.py --annot annotations.txt ../IQ-TREE2/final_aa_trimmed_per_gene.aln.treefile mags_phylogeny.xml
graphlan.py mags_phylogeny.xml mags_phylogeny.svg --dpi 200 --size 13 --pad 0.01 --external_legends
```



