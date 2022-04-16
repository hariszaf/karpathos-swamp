#!/usr/bin/env python

import sys 

aa_file = open(sys.argv[1])

bin_seq = {}
for line in aa_file: 

	line        = line.split("\t")
	bin_number  = line[0]
	seq         = line[1][:-1]
	bin_seq[bin_number] = seq
	seq_len = len(seq)

ofile = sys.argv[1] + ".aln"
#ofile = "numbers.tsv"
o     = open(ofile, "w")

for i in range(1,290):

	bin_name = "bin" + str(i)	
	try: 
		o.write(bin_seq[bin_name] + "\n")
	except:
		o.write("-"*seq_len + "\n")
#	o.write(">bin_" + str(i) + "\n")
