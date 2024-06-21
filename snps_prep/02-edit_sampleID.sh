#!/bin/bash

# Script for editing labels in the vcf, shortening sample ID to remove unncessary characters/paths in sample IDs

# Example: I want to shorten sample ID from "/edegreef/whales/BOW_bam/88_Pang_S51" to "88_Pang"
# Can also go straight to "bcftools reheader" if want to use a pre-made sample ID list.

variants=bowhead_RWmap_allvariants

#1) Change sample ID in vcf to just contain the individual ID (removing the path that is currently in sample name).

# Zip and make index file if not done already
bgzip -c $variants.vcf > $variants.vcf.gz
tabix -p vcf $variants.vcf.gz

# Make list of current sample ID names
vcf-query -l $variants.vcf.gz > vcf_sampleID

# Shorten sample IDs by removing text up to the last "/"
awk '{print $NF}' FS=/ vcf_sampleID > vcf_sampleID_shortened

# also want to remove last 4 characters (e.g. _S43 or something from sequencing lane), not related to sample ID
sed 's/....$//' vcf_sampleID_shortened > vcf_sampleID_shortened2

# Modify header in vcf file to change sample ID to the shortened ID list
/home/degreefe/programs/bcftools-1.9/bcftools reheader -s vcf_sampleID2 $variants.vcf.gz -o $variants.ID.vcf.gz
tabix -p vcf $variants.ID.vcf.gz

