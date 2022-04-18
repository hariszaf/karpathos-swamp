#!/usr/bin/env python

import os, sys


data = open("data.csv", "r")
annotation_file = open("annotations.txt", "w")

annotation_file.write(
"title	MAGs phylogeny\n\
title_font_size	25\n\
total_plotted_degrees	340\n\
start_rotation	270\n\
class_legend_font_size	12\n\
annotation_legend_font_size	11\n\
clade_marker_size	3.0\n\
clade_separation	0.35\n\
annotation_background_separation	0.1\n\
annotation_background_offset	0.05\n\
annotation_background_width	0.1\n\
ignore_branch_len	0\n"
)


annotation_file.write("ring_internal_separator_thickness	1	0.5\n\
ring_separator_color	1	#696969\n\
ring_label	1	Phylum\n\
ring_label_color	1	#696969\n")

# EMPTY RING
annotation_file.write("ring_internal_separator_thickness	2	0.5\n\
ring_separator_color	2	#696969\n\
ring_label_color	2	#FFFFFF\n")


# SAMPLES' ABUNDANCES
annotation_file.write("ring_internal_separator_thickness	3	0.5\n\
ring_separator_color	3	#696969\n\
ring_label	3	top sediment\n\
ring_label_color	3	#696969\n")

annotation_file.write("ring_internal_separator_thickness	4	0.5\n\
ring_separator_color	4	#696969\n\
ring_label	4	bottom sediment\n\
ring_label_color	4	#696969\n")

annotation_file.write("ring_internal_separator_thickness	5	0.5\n\
ring_separator_color	5	#696969\n\
ring_label	5	orange aggreg\n\
ring_label_color	5	#696969\n")

annotation_file.write("ring_internal_separator_thickness	6	0.5\n\
ring_separator_color	6	#696969\n\
ring_label	6	pink aggreg\n\
ring_label_color	6	#696969\n")

annotation_file.write("ring_internal_separator_thickness	7	0.5\n\
ring_separator_color	7	#696969\n\
ring_label	7	combined\n\
ring_label_color	7	#696969\n")

annotation_file.write("ring_internal_separator_thickness	8	0.5\n\
ring_separator_color	8	#696969\n\
ring_label	8	orange carpet\n\
ring_label_color	8	#696969\n")



colors = [
        "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
        "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
        "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
        "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD",
        "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
        "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
        "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
        "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
        "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
        "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
        ]


phyla_colors = {}
phyla_counter = 0 
max_abundance = 0
for line in data:

    line = line.split(",")

    bin_name    = line[0]
    
    abundance_1 = line[1]
    abundance_2 = line[2]
    abundance_3 = line[3]
    abundance_4 = line[4]
    abundance_5 = line[5]
    abundance_6 = line[6]

    phylum      = line[8]
    complet     = line[14]

    if phylum not in phyla_colors: 
        phylum_color = colors[phyla_counter]
        phyla_colors[phylum] = phylum_color
        phyla_counter += 1 

    # TAXONOMY
    annotation_file.write(bin_name + "\t" + "ring_color\t" + "1" + "\t" + phyla_colors[phylum] + "\n")
    annotation_file.write(bin_name + "\t" + "ring_width\t" + "1" + "\t" + "1.0" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height\t" + "1" + "\t" + "0.3" + "\n")    

    # LEAVE EMPTY RING
    annotation_file.write(bin_name + "\t" + "ring_color\t" + "2" + "\t" + "#FFFFFF" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_width\t" + "2" + "\t" + "1.0" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height\t" + "2" + "\t" + "0.3" + "\n")    

    # ABUNDANCE - SAMPLE 1
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "3" + "\t" + str(abundance_1) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "3" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "3" + "\t" + "0.3" + "\n")

    # ABUNDANCE - SAMPLE 2
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "4" + "\t" + str(abundance_2) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "4" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "4" + "\t" + "0.3" + "\n")

    # ABUNDANCE - SAMPLE 3
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "5" + "\t" + str(abundance_3) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "5" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "5" + "\t" + "0.3" + "\n")

    # ABUNDANCE - SAMPLE 4
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "6" + "\t" + str(abundance_4) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "6" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "6" + "\t" + "0.3" + "\n")

    # ABUNDANCE - SAMPLE 5
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "7" + "\t" + str(abundance_5) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "7" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "7" + "\t" + "0.3" + "\n")

    # ABUNDANCE - SAMPLE 6
    annotation_file.write(bin_name + "\t" + "ring_alpha" + "\t" + "8" + "\t" + str(abundance_6) + "\n")
    annotation_file.write(bin_name + "\t" + "ring_color" + "\t" + "8" + "\t" + "#AA00AA" + "\n")
    annotation_file.write(bin_name + "\t" + "ring_height"+ "\t" + "8" + "\t" + "0.3" + "\n")

    # BIN NAME 
    annotation_file.write(bin_name + "\t" + "annotation" + "\t" + bin_name + "\n")
    annotation_file.write(bin_name + "\t" + "annotation_rotation" + "\t" + "90" + "\n")
    annotation_file.write(bin_name + "\t" + "annotation_background_color" + "\t" + "#f3f6f4" + "\n")
    annotation_file.write(bin_name + "\t" + "annotation_background_edge_color" + "\t" + "#f3f6f4" + "\n")

    # CHECK COMPLET
    if "1" in complet:
        annotation_file.write(bin_name + "\t" "clade_marker_shape" + "\t" + "*" + "\n")
        annotation_file.write(bin_name + "\t" "clade_marker_size" + "\t" + "30" + "\n")
        annotation_file.write(bin_name + "\t" "clade_marker_color" + "\t" + "#2c8f4a" + "\n")

    # CHECK FOR NOTABLE TAXA
    if phylum == "p__Afoulota":

        annotation_file.write(bin_name + "\t" "annotation" + "\t" + bin_name + "|" + phylum  + "\n")
        annotation_file.write(bin_name + "\t" + "annotation_rotation" + "\t" + "90" + "\n")
        annotation_file.write(bin_name + "\t" + "annotation_background_color" + "\t" + "grey" + "\n")
        annotation_file.write(bin_name + "\t" + "annotation_background_edge_color" + "\t" + "grey" + "\n")


    # MAX ABUNDANCE VALUE 
    max_bin_abund = max(float(abundance_1), float(abundance_2), float(abundance_3), float(abundance_4), float(abundance_5), float(abundance_6))
    if max_bin_abund > max_abundance: 
        max_abundance = max_bin_abund


for phylum, color in phyla_colors.items(): 

    annotation_file.write(phylum + "\t" + "annotation" + "\t" + phylum + "\n")
    annotation_file.write(phylum + "\tclade_marker_color\t" + color + "\n")
    annotation_file.write(phylum + "\tclade_marker_size\t" + "75\n")


annotation_file.write(">90% completeness" + "\tclade_marker_shape\t" + "*" + "\n")
annotation_file.write(">90% completeness" + "\tclade_marker_size\t" + "75\n")
annotation_file.write(">90% completeness" + "\tclade_marker_color\t" + "#2c8f4a\n")


# ABUNDANCE RANGE
step = max_abundance / 10
palette = ["#FFFFFF", "#F5E2F5", "#ECC6EC", "#E2AAE2", "#D98DD9", "#CF71CF", "#C655C6", "#BC38BC", "#B31CB3", "#AA00AA" ]

for i in range(10): 
    annotation_file.write(str(int(i*step)) + "\tannotation\t" + "abund:" + str(int(i*step)) + "\n")
    annotation_file.write(str(int(i*step)) + "\tclade_marker_shape\t" + "s" + "\n")
    annotation_file.write(str(int(i*step)) + "\tclade_marker_size\t" + "75\n")
    annotation_file.write(str(int(i*step)) + "\tclade_marker_color\t" + palette[i] + "\n")


bs_clades = open("clade_names.txt", "r")
for clade in bs_clades: 

    uf_bs_clade = clade.split("/")[1]
    uf_bs_value = float(uf_bs_clade)

    if uf_bs_value <70:
        annotation_file.write(clade[:-1] + "\t" + "clade_marker_color" + "\t" + "b" + "\n")

    elif uf_bs_value < 90: 
        annotation_file.write(clade[:-1] + "\t" + "clade_marker_color" + "\t" + "g" + "\n")

    else: 
        annotation_file.write(clade[:-1] + "\t" + "clade_marker_color" + "\t" + "k" + "\n")



annotation_file.write("bootstrap < 70" + "\tannotation\t" + "bootstrap < 70" + "\n")
annotation_file.write("bootstrap < 70" + "\tclade_marker_shape\t" + "_" + "\n")
annotation_file.write("bootstrap < 70" + "\tclade_marker_size\t" + "105\n")
annotation_file.write("bootstrap < 70" + "\tclade_marker_color\t" + "b" + "\n")

annotation_file.write("bootstrap < 90" + "\tannotation\t" + "bootstrap < 90" + "\n")
annotation_file.write("bootstrap < 90" + "\tclade_marker_shape\t" + "_" + "\n")
annotation_file.write("bootstrap < 90" + "\tclade_marker_size\t" + "105\n")
annotation_file.write("bootstrap < 90" + "\tclade_marker_color\t" + "g" + "\n")


annotation_file.write("bootstrap > 90" + "\tannotation\t" + "bootstrap > 90" + "\n")
annotation_file.write("bootstrap > 90" + "\tclade_marker_shape\t" + "_" + "\n")
annotation_file.write("bootstrap > 90" + "\tclade_marker_size\t" + "105\n")
annotation_file.write("bootstrap > 90" + "\tclade_marker_color\t" + "k" + "\n")




# #FFFFFF --> (255, 255, 255)
# #F5E2F5 --> (245, 226, 245)
# #ECC6EC --> (236, 198, 236)
# #E2AAE2 --> (226, 170, 226)
# #D98DD9 --> (217, 141, 217)
# #CF71CF --> (207, 113, 207)
# #C655C6 --> (198, 85, 198)
# #BC38BC --> (188, 56, 188)
# #B31CB3 --> (179, 28, 179)
# #AA00AA --> (170, 0, 170)

#palette = [(255, 255, 255), (245, 226, 245), (236, 198, 236), (226, 170, 226), (217, 141, 217), (207, 113, 207), (198, 85, 198), (188, 56, 188), (179, 28, 179), (170, 0, 170)]
