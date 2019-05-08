#!/usr/bin/env python
'''This script takes a VCF file and counts the numbers of heterozygous and missing genotypes per site. Written by Verena Kutschera, January 2017.

Usage:
python countHeterozygousGenotypes.py myData.vcf.gz het_missing_sites.out'''

import sys
import gzip
import re

vcfRead = gzip.open(sys.argv[1], 'r') #open the gzipped VCF file.
output_file = open(sys.argv[2], 'w') #write to outputfile (unzipped).
output_file.write('scaffold\tposition\tno_het_genotypes\tprop_het_genotypes\tno_missing_genotypes\n')

lines = vcfRead.readlines()

for site in lines: #loop through sites in vcf file. For each site count the number of missing and heterozygous sites.
    if not site.startswith('#'):
        col = site.strip().split()
        scaffold = str(col[0])
        position = int(col[1])
        het = []
        missing = []
        for j in col:
            geno = j.split(':')[0]
            if geno == './.':
                missing.append(1)
            if re.match("0/[1-9]", geno, flags=0):
                het.append(1)
            if re.match("[1-9]/0", geno, flags=0):
                het.append(1)
            if re.match("1/[2-9]", geno, flags=0):
                het.append(1)
            if re.match("2/[3-9]", geno, flags=0):
                het.append(1)
            if re.match("3/[4-9]", geno, flags=0):
                het.append(1)
            if re.match("4/[5-9]", geno, flags=0):
                het.append(1)
            if re.match("5/[6-9]'", geno, flags=0):
                het.append(1)
            if re.match("6/[7-9]'", geno, flags=0):
                het.append(1)
        n_het = sum(het)
        n_missing = sum(missing)
        length = len(col)
        n_inds = float(length) - 9.0 - float(n_missing)
        prop_het = float(n_het) / float(n_inds)
        output_file.write('%s\t%i\t%i\t%f\t%i\n' % (scaffold, position, n_het, prop_het, n_missing))
            
vcfRead.close()
output_file.close()
