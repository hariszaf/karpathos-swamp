#!/bin/bash

for bin in `awk '{print $1}' fa.list`
do

    echo $bin

    grep ">" $bin | sed 's/>//' > tmp

    awk -v var=$bin 'BEGIN{ ans=var} { print $0 "\t" ans}' tmp >> external_binning_of_contigs.txt

done


