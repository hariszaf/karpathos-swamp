#! /usr/bin/gawk -f

###############################################################################
#
# usage: ./mapConts.awk map_int_cont_name_to_anvios.tsv 
#                       external_binning_of_contigs.txt > external_binning_of_contigs_mapped.txt 
#
###############################################################################
BEGIN {
    FS="\t"
}
# Load the data in associative arrays.

(ARGIND==1) {

    #initiate an array with the desired NCBI ids to include only microbes.
    init_contig_to_anvio[$2]=$1
}
#Load all the database pairs files from all sources and channels of PREGO
(ARGIND==2){


    bin[$1]=$2

}
END{

    for (entry in bin){
        print init_contig_to_anvio[entry] "\t" bin[entry]
    }
}

