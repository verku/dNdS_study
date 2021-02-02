# dNdS_study
Custom python scripts used in bioinformatics analysis for Mugal et al., Molecular Biology and Evolution, Volume 37, Issue 1, January 2020, Pages 260–279, https://doi.org/10.1093/molbev/msz203


# bootstrap_codon_alignment.py
- Creates one bootstrap replicate of a multiple alignment of protein-coding genes in FASTA format.
- Randomly draws triplets of bases (codons) with replacement, until the original length of the alignment is reached.
- The bootstrap replicate is written to a separate FASTA file, into a directory with the same name as the original alignment.
- Requires python3 and biopython 1.72. Not optimized for large alignments.

# fasta_N_nucleotide_content.py
- Counts the numbers of the different nucleotides incl. Ns, as well as the proportion of Ns, of polymorphic sites coded as IUPAC codes and the GC content (only taking Gs and Cs into account) in a FASTA file.
- Requires python2.7 and  Biopython 1.43 onwards.

# countHeterozygousGenotypes.py
- Takes a gzipped VCF file and counts the numbers of heterozygous and missing genotypes per site.
- The output is printed to a table.
- Requires python2.7.

# findCpGpolymorphism.py
- Takes a gzipped VCF file containing variant and monomorphic/invariant sites, and finds the positions of CpG, CpA and TpG sites, incl CpG/CpA and TpG/CpG polymorphisms within the resequencing data.
- Not only REF and ALT, but also allele frequency (AF) and genotype information (individual GT) are considered. 
- Site pairs with positions not directly adjacent to each other will be ignored. 
- The output is printed to a bed file.
- Written for VCF files generated with GATK 3.4.0.
- Requires python2.7.
