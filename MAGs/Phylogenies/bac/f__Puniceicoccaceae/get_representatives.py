#!/usr/bin/env python3

import glob
import os
import subprocess
import sys
from subprocess import call
import HTSeq

dir_path = os.getcwd() + '/'

representatives = open("repr_accessions_numbers.tsv", "r").readlines()


genome_file = "mags/bin63-contigs.fa"
bashCommand = 'prokka --centre C --locustag L --cpus 20 '+ genome_file+' -o ' + dir_path + '/PROKKA-annotation --prefix Isolate_' + genome_file.split("/")[1]
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


# BLAST QUERY MAG AGAINST CLOSE RELATIVES
genomes = []
for genome in glob.glob(dir_path + '/CLOSE_RELS/*.faa'):
    genomes.append(genome)


aimed_for = ''
for genome in glob.glob(dir_path + '/PROKKA-annotation/*.faa'):
    genomes.append(genome)
    aimed_for = genome

raw_data = {}
for FILE1 in genomes:
    temp = {}
    for FILE2 in genomes:
        if aimed_for in [FILE1, FILE2]:
            bashCommand = 'blastp -num_threads 20 -evalue 0.00001 -qcov_hsp_perc 50.0 -subject ' + FILE2 + ' -query ' + FILE1 + ' -outfmt 6  > ' + dir_path  + '/temp.m8'
            notneeded = call(bashCommand, shell=True)
            total = 0
            for read in HTSeq.FastaReader(FILE1):
                total +=1
            matches = []
            for line in open(dir_path  + '/temp.m8'):
                timber = line.split('\t')
                if float(line.split('\t')[2]) > 40.0:
                    if line.split('\t')[0] in matches:
                        continue
                    else:
                        matches.append(line.split('\t')[0])
            temp[FILE2] = [total,len(matches)]
    raw_data[FILE1] = temp

# NOW COUNT AND WRITE
number_tested = 0
same_genus = 0
different_genus = 0

outputting_pocp = open(dir_path + '/POCP_results.tab','w')
outputting_pocp.write('#Matching species\tPOCP (%)\n')
assigned_genus = ''
assigned_pocp = 0.0
for genome1 in glob.glob(dir_path + '/PROKKA-annotation/*.faa'):
    for genome2 in genomes:
        if genome1 != genome2:
            g1d = raw_data[genome1][genome2]
            g2d = raw_data[genome2][genome1]
            C1 = g1d[1]
            C2 = g2d[1]
            T1 = g1d[0]
            T2 = g2d[0]
            PCOP = ((C1+float(C2))/(T1+float(T2)))*100
            number_tested +=1
            if PCOP > 50.0:
                same_genus +=1
                outputting_pocp.write('POCP suggests the input genome and '+ genome2.split('/')[-1:][0].split('_genomic')[0] + ' belong to the same genus with a POCP was;\t' + str(PCOP) + '%\n')
                assigned_genus = genome2.split('/')[-1:][0].split('_genomic')[0].split('_')[0]
                assigned_pocp = PCOP
            else:
                outputting_pocp.write('POCP suggests seperate genus between the input genome and '+ genome2.split('/')[-1:][0].split('_genomic')[0] + ' as POCP was;\t' + str(PCOP) + '%\n')


#output.write('The input genome was assigned to ' + assigned_genus +' with a POCP value of ' + str(assigned_pocp) +'%\n')
