#!/bin/bash

# Script for editing sample IDs in the vcf 
# Example: I want to shorten sample ID from "/edegreef/whales/BOW_bam/88_Pang_S51" to "88_Pang"
# Can also go straight to "bcftools reheader" if want to use a pre-made sample ID list.

snps=BOW_allvariantcalls_pp

# Zip and index vcf if not done already
bgzip -c $snps.vcf > $snps.vcf.gz
tabix -p vcf $snps.vcf.gz

# Make list of current sample ID names
vcf-query -l $snps.vcf.gz > vcf_sampleID

# Shorten sample IDs by removing text up to the last "/"
awk '{print $NF}' FS=/ vcf_sampleID > vcf_sampleID_shortened

# I also want to remove last 4 characters in my sample ID (e.g. _S43 or something from sequencing lane), not related to sample ID
sed 's/....$//' vcf_sampleID_shortened > vcf_sampleID_shortened2

# Modify header in vcf file to change sample ID to the shortened ID list
bcftools reheader -s vcf_sampleID_shortened2 $snps.vcf.gz -o $snps.ID.vcf.gz

# Index new vcf file
tabix -p vcf $snps.ID.vcf.gz
