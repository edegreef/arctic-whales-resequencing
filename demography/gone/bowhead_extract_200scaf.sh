#!/bin/bash


# Set up these before running
snps_prefix=bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.subsamp11
bcftools_path=/home/degreefe/programs/bcftools-1.9
scaf_list=top200_bowheadscaf_list.txt #list of scaffolds to keep


## want to keep just the scaffolds in the list

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

bgzip $snps_prefix.vcf
tabix -p vcf $snps_prefix.vcf.gz

# Filter vcf for these scaffolds
#$bcftools_path/bcftools filter --regions $list $snps_prefix.vcf.gz > $snps_prefix.top200scaf.vcf
$bcftools_path/bcftools filter --regions $list $snps_prefix.vcf.gz > $snps_prefix.top200scaf.vcf

# Can make a list of scaffolds to double check
#grep -v "^#" $snps_prefix.top200scaf.vcf | cut -f1 | sort | uniq > checktop200.txt

