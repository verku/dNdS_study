# dNdS_study
Custom python scripts used in bioinformatics analysis for Mugal et al.

# bootstrap_codon_alignment.py
- Creates one bootstrap replicate of a multiple alignment of protein-coding genes in FASTA format.
- Randomly draws triplets of bases (codons) with replacement, until the original length of the alignment is reached.
- The bootstrap replicate is written to a separate FASTA file, into a directory with the same name as the original alignment.
- Requires python3 and biopython 1.72. Not optimized for large alignments.

# fasta_N_nucleotide_content.py
- Counts the numbers of the different nucleotides incl. Ns, as well as the proportion of Ns, of polymorphic sites coded as IUPAC codes and the GC content (only taking Gs and Cs into account). 
- Requires python2.7 and  Biopython 1.43 onwards.
