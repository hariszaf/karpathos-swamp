#!/bin/bash


for bin  in `awk '{print $1}' bac.list`
do

barrnap -k bac --outseq rrnas_$bin < $bin  > bac_rrnas/rrnas_$bin.gff

done


#for bin  in `awk '{print $1}' arc.list`
#do

#barrnap -k bac --outseq rrnas_$bin.fa < $bin  > arc_rrnas/rrnas_$bin.gff

#done
