#!/bin/bash

# Script to make a vcf containing select scaffolds with a list. 
# I am filtering out small scaffolds, using a 100 kb as the threshold. List created by looking at reference genome scaffold lengths.

# Input snp prefix, and name of file that contains a list of scaffolds to include
snps=BOW_SNPS.filter1.miss.biallel
scaf_list=scaf_min100kb

# Index vcf if not indexed already
bgzip $snps.vcf
tabix -p vcf $snps.vcf.gz

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

# Filter vcf for these scaffolds
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list $snps.vcf.gz > $snps.min100kb.vcf

# Can make a list of scaffolds to double check
grep -v "^#" $snps.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt