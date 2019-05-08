#!/usr/bin/env python3

'''
Python script to create one bootstrap replicate of a multiple alignment of protein-coding genes in FASTA format.
Randomly draws triplets of bases (codons) with replacement, until the original length of the alignment is reached.
The bootstrap replicate is written to a separate FASTA file, into a directory with the same name as the original alignment.

Requires python3 and biopython 1.72. Not optimized for large alignments.

Usage:
./bootstrap_codon_alignment.py [path to original alignment in fasta format] [path to output alignment in fasta format]
'''

import sys
import os
from Bio import AlignIO
from Bio import SeqIO
import Bio.Align
import random


input_file = open(sys.argv[1], 'r') # take input file name from command line
output_file = open(sys.argv[2], 'w') # path to output file

def triplet_idx_list(fasta):
	aln = AlignIO.read(input_file, "fasta")
	length = int(aln.get_alignment_length()) # length of the alignment
	idxlist = [] # list to store every third index
	if length % 3 == 0: # only continue if the alignment is a multiple of three
		for i in range(0,length): # up to the length of the alignment
			if i % 3 == 0: # append index to list if division by three possible (incl. zero)
				idxlist.append(i)
	return idxlist, aln

def random_aln_generator(idxlist, aln):
	randomlist = random.choices(idxlist, k=len(idxlist)) # draw randomly from the index list of codon start positions with replacement
	randomalign = aln[:, randomlist[0]:randomlist[0]+3] # list with random triplets, each stored as alignment object
	for i in randomlist[1:]: # loop through random list of indices of codon start positions 
		j = i+3 # create range for each random triplet
		randomalign += aln[:, i:j]
	return randomalign


# create desired number of bootstrap replicates in fasta format
indexlist, alignment = triplet_idx_list(input_file) # get the indices for the original alignment 

replicate = random_aln_generator(indexlist, alignment)
SeqIO.write(replicate, output_file, "fasta")

input_file.close()
output_file.close()