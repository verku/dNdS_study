#!/usr/bin/env python
'''This script takes a gzipped VCF file containing variant and monomorphic/invariant sites, 
and finds the positions of CpG, CpA and TpG sites, incl CpG/CpA and TpG/CpG polymorphisms within the resequencing data.
Not only REF and ALT, but also allele frequency (AF) and genotype information (individual GT) are considered.
Site pairs with positions not directly adjacent to each other will be ignored. The output is printed to a bed file.
Written for python 2.7, by Verena Kutschera, May 2017.

Specifically, the script: 
--> Finds sites with "C" fixed
--> if yes, checks if the allele at the following site is "G" and/or "A" (positions must be subsequent, no gaps allowed)
--> if yes, prints the positions of the two sites into a new bed file

--> Finds sites with "C" and/or "T"
--> if yes, checks if allele at the following site is fixed for "G" (positions must be subsequent, no gaps allowed)
--> if yes, prints the positions of the two sites into a new bed file

Usage:
python findCpGpolymorphism.py resequencingData.vcf.gz CpG_polymorphism.bed'''

import gzip
import sys
from collections import deque

vcfRead = gzip.open(sys.argv[1], 'r') #open a gzipped VCF file
bedWrite = open(sys.argv[2], 'w') #write to outputfile in bed format (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

### a set of functions to be called in the for loop through the VCF file
# splits each line and returns a list of arguments taken by the other functions
def getSiteInfo(site):
	col = site.strip().split()
	chrom = col[0] # the scaffold name
	pos = col[1] # the position
	ref = col[3] # the reference base
	alt = col[4] # the alt base
	info = col[7].strip().split(";") # all info for the site in one list
	geno = [] # all genotypes in a list to be checked by the functions
	for gt in col[9:]:
		geno.append(gt.split(':')[0])
	return chrom, pos, ref, alt, info, geno

# finds sites fixed for "C" or "G", and return scaffold name, position and fixed allele
def CorGfixed(chrom, pos, ref, alt, info, geno):
	if ref=="C" or ref=="G":
		if alt==".": # definitely monomorphic for ref allele
			return chrom, pos, ref # "C_or_G_fixed_ref"
		elif (geno.count("0/0")+geno.count("./.")==len(geno)) and any("AF=0.00" in e for e in info): # any allele as alt, but all genotypes "0/0"
			return chrom, pos, ref # "C_or_G_fixed_alt_af0"
		elif (alt=="C" or alt=="G") and (geno.count("1/1")+geno.count("./.")==len(geno)) and any("AF=1.00" in e for e in info): # alt fixed for "C" or "G"
			return chrom, pos, alt # "C_or_G_fixed_alt_af1"
	elif alt=="C" or alt=="G": # capturing sites that are monomorphic for the alt allele but have a different allele as ref
		if geno.count("1/1")+geno.count("./.")==len(geno) and any("AF=1.00" in e for e in info): # site monomorphic for alt allele
			return chrom, pos, alt # "C_or_G_fixed_alt_af1"

# finds sites with "G" and/or "A" (fixed or polymorphic)
def GApol(chrom, pos, ref, alt, info, geno):
	if ref=="A" or ref=="G": # sites fixed for "A" or "G"
		if alt==".": # definitely monomorphic for ref allele
			return chrom, pos, ref # "G_or_A_fixed_ref"
		elif (geno.count("0/0")+geno.count("./.")==len(geno)) and any("AF=0.00" in e for e in info): # any allele as alt, but all genotypes "0/0"
			return chrom, pos, ref # "G_or_A_fixed_af0"
		elif alt=="A" or alt=="G" or alt=="G,A" or alt=="A,G":
			return chrom, pos, ref + "," + alt # "G_or_A_ref_alt"
	elif alt=="A" or alt=="G" or alt=="G,A" or alt=="A,G":
		if ref!="A" or ref!="G": # both conditions below target the same case: if the ref allele is not "A" or "G", it is not allowed among the genotypes
			if geno.count("1/1")+geno.count("1/2")+geno.count("2/2")+geno.count("2/1")+geno.count("./.")==len(geno): # to make sure only the alt alleles are represented among the genotypes
				return chrom, pos, alt # "G_or_A_alt"
			elif any("0" in g for g in geno): # no "0" allowed in any of the genotypes
				return None
		elif ref=="A" or ref=="G": # get "A"/"G" polymorphism
			return chrom, pos, ref + "," + alt # "G_or_A_alt_ref"

# finds sites with "C" and/or "T" (fixed or polymorphic)
def CTpol(chrom, pos, ref, alt, info, geno):
	if ref=="T" or ref=="C":
		if alt==".": # definitely monomorphic for ref allele
			return chrom, pos, ref # "C_or_T_fixed_ref"
		elif (geno.count("0/0")+geno.count("./.")==len(geno)) and any("AF=0.00" in e for e in info): # any allele as alt, but all genotypes "0/0"
			return chrom, pos, ref # "C_or_T_fixed_af0"
		elif alt=="T" or alt=="C" or alt=="C,T" or alt=="T,C":
			return chrom, pos, ref + "," + alt # "C_or_T_ref_alt"
	elif alt=="T" or alt=="C" or alt=="C,T" or alt=="T,C":
		if ref!="T" or ref!="C": # both conditions below target the same case: if the ref allele is not "T" or "C", it is not allowed among the genotypes
			if geno.count("1/1")+geno.count("1/2")+geno.count("2/2")+geno.count("2/1")+geno.count("./.")==len(geno): # to make sure only the alt alleles are represented among the genotypes
				return chrom, pos, alt # "C_or_T_alt"
			elif any("0" in g for g in geno): # no "0" allowed in any of the genotypes
				return None
		elif ref=="T" or ref=="C": # should not happen because of first if clause, but captures any "A"/"G" polymorphism but no other allele present at that site
			return chrom, pos, ref + "," + alt # "C_or_T_alt_ref"

# converts line from VCF into bed format
def bedFormat(chrom, pos, base):
	startpos = int(pos)-1
	return str(chrom) + "\t" + str(startpos) + "\t" + str(pos) + "\t" + str(base)
### end of functions


# Loop through VCF file, store subsequent lines in a deque to be checked by the functions
lineDeque = deque(maxlen=2) # store the deque of a maximum length of two
for line in vcfRead:
	if line.startswith('#'): # skip the header section
		continue
	else:
		lineDeque.append(line) # move the sliding window one line forward
		if len(lineDeque)==2: # only continue if there are exactly two lines in the deque
			siteInfo = getSiteInfo(lineDeque[0]) # split the first line in the deque into its parts to get the relevant info needed for the other functions
			fixedSite = CorGfixed(*siteInfo) # run the function to test if "C" or "G" are fixed at that site. The arguments are taken from the list siteInfo, indicted by a star
			CTpolSite = CTpol(*siteInfo) # run the function to test if there is a "C" or "T" at that site. The arguments are taken from the list siteInfo, indicted by a star

			# find CpG or CpA sites in the resequencing data, by checking the positions of the two respective alleles
			if fixedSite!=None and fixedSite[2]=="C": # if the site is fixed for "C"
				nextSiteInfo = getSiteInfo(lineDeque[1]) # split the next site so that it's parameters can be taken by the functions
				if nextSiteInfo[0]==siteInfo[0] and (int(nextSiteInfo[1])-1==int(siteInfo[1])): # check that the sites are on the same scaffold and consecutive
					GApolNextSite = GApol(*nextSiteInfo) # run the function to check if the next site is a "G", "A", "G,A", or "A,G"
					if GApolNextSite!=None: # check if the next site is a "G", "A", "G,A", or "A,G"
						print >> bedWrite, bedFormat(*fixedSite), "\n", bedFormat(*GApolNextSite) # print the sites in bed format

			# find CpG or TpG sites in the resequencing data, by checking the positions of the two respective alleles
			elif CTpolSite!=None: # if alt allele is "C", "T", "C,T" or "T,C"
				nextSiteInfo = getSiteInfo(lineDeque[1]) # split the next site so that it's parameters can be taken by the functions
				if nextSiteInfo[0]==siteInfo[0] and (int(nextSiteInfo[1])-1==int(siteInfo[1])): # check that the sites are on the same scaffold and consecutive
					fixedNextSite = CorGfixed(*nextSiteInfo) # run the function to check if the next site is fixed for "C" or "G"
					if fixedNextSite!=None and fixedNextSite[2]=="G": # check that the next site is fixed for "G"
						print >> bedWrite, bedFormat(*CTpolSite), "\n", bedFormat(*fixedNextSite) # print the sites in bed format

vcfRead.close()
bedWrite.close()
