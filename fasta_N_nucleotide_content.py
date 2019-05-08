#!/usr/bin/env python2

'''This script counts the numbers of the different nucleotides incl. Ns, as well as the proportion of Ns, of polymorphic sites coded as IUPAC codes and the GC content (only taking Gs and Cs into account). 
It requires Biopython 1.43 onwards.

Usage:
python fasta_N_nucleotide_content.py input.fa output.table

'''

import sys

input_file = open(sys.argv[1], 'r') #take input file name from command line
output_file = open(sys.argv[2],'w') #take output file name from command line
output_file.write('Gene\tLength\tN\tN%\tpolymorphic%\tA\tC\tG\tT\tGC%\n')

from Bio import SeqIO

for cur_record in SeqIO.parse(input_file, "fasta") :
	records = cur_record.upper()
	gene_name = cur_record.name
	N_count = records.seq.count('N')
	A_count = records.seq.count('A')
	C_count = records.seq.count('C')
	G_count = records.seq.count('G')
	T_count = records.seq.count('T')
	R_count = records.seq.count('R')
	Y_count = records.seq.count('Y')
	S_count = records.seq.count('S')
	W_count = records.seq.count('W')
	K_count = records.seq.count('K')
	M_count = records.seq.count('M')
	B_count = records.seq.count('B')
	D_count = records.seq.count('D')
	H_count = records.seq.count('H')
	V_count = records.seq.count('V')
	length = len(records.seq)
	n_percentage = float(N_count) / length
	poly_percentage = float(R_count + Y_count + S_count + W_count + K_count + M_count + B_count + D_count + H_count + V_count) / length
	gc_percentage = float(C_count + G_count) / length
	output_line = '%s\t%i\t%i\t%f\t%f\t%i\t%i\t%i\t%i\t%f\n' % \
	(gene_name, length, N_count, n_percentage, poly_percentage, A_count, C_count, G_count, T_count, gc_percentage)
	output_file.write(output_line)

output_file.close()
input_file.close() 
