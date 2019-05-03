This repository contains python scripts to handle FASTA files.

All scripts were written in python3.

# bootstrap_codon_alignment.py
- Creates one bootstrap replicate of a multiple alignment of protein-coding genes in FASTA format.
- Randomly draws triplets of bases (codons) with replacement, until the original length of the alignment is reached.
- The bootstrap replicate is written to a separate FASTA file, into a directory with the same name as the original alignment.
- Requires python3 and biopython 1.72. Not optimized for large alignments.

